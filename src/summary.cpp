//
//  summary.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/17.
//  Copyright © 2024 Schloss Lab. All rights reserved.
//

#include "summary.h"
#include "sequence.h"

//**************************************************************************
Summary::Summary(int p) {
    processors = p;
    total = 0;
    numUniques = 0;
    hasCount = false;
}
//**************************************************************************
vector<long long> Summary::getDefaults() {
        vector<long long> locations;

        //number of sequences at 2.5%
        long long ptile0_25	= 1+(long long)(total * 0.025);
        //number of sequences at 25%
        long long ptile25		= 1+(long long)(total * 0.250);
        long long ptile50		= 1+(long long)(total * 0.500);
        long long ptile75		= 1+(long long)(total * 0.750);
        long long ptile97_5	= 1+(long long)(total * 0.975);
        long long ptile100	= (long long)(total);

        locations.push_back(1);
        locations.push_back(ptile0_25);
        locations.push_back(ptile25);
        locations.push_back(ptile50);
        locations.push_back(ptile75);
        locations.push_back(ptile97_5);
        locations.push_back(ptile100);

        return locations;
}
//**************************************************************************
long long Summary::getValue(map<int, long long>& spots, double value) {

    long long percentage = 1+(long long)(total * value * 0.01);
    long long result = 0;
    long long totalSoFar = 0;
    long long lastValue = 0;

    //minimum
    if ((spots.begin())->first == -1) { result = 0; }
    else {result = (spots.begin())->first; }

    for (auto it = spots.begin(); it != spots.end(); it++) {
        long long value = it->first; if (value == -1) { value = 0; }
        totalSoFar += it->second;

        if (((totalSoFar <= percentage) && (totalSoFar > 1)) ||
            ((lastValue < percentage) && (totalSoFar > percentage))){
            result = value;
        }
        lastValue = totalSoFar;
    }

    return result;

}
//**************************************************************************
vector<int> Summary::getValues(map<int, long long>& positions) {
        vector<long long> defaults = getDefaults();
        vector<int> results; results.resize(7,0);
        long long meanPosition; meanPosition = 0;
        long long totalSoFar = 0;
        int lastValue = 0;

        // minimum
        if ((positions.begin())->first != -1) {
            results[0] = (positions.begin())->first;
        }

        results[1] = results[0]; results[2] = results[0];
        results[3] = results[0]; results[4] = results[0];
        results[5] = results[0];

        for (auto it = positions.begin(); it != positions.end(); it++) {

            int value = it->first; if (value == -1) { value = 0; }
            meanPosition += (value*it->second);
            totalSoFar += it->second;

            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) ||
                ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){
                results[1] = value;
            } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||
                ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) {
                results[2] = value;
            } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||
                ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {
                results[3] = value;
            } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||
                ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {
                results[4] = value;
            } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||
                ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {
                results[5] = value;
            } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {
                results[6] = value;
            } //save value
            lastValue = totalSoFar;
        }
        // maximum
        results[6] = (positions.rbegin())->first;


        double meansPosition = meanPosition / (double) total;

        // mean
        results.push_back(meansPosition);

        return results;
}
//**************************************************************************
struct seqSumData {
    // fasta
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    map<int, long long> numNs;

    vector<int> counts; // abundances for each read, all ones if no count info
    vector<vector<int>> summary;

    double start;
    double end;

    seqSumData(){
        start = 0;
        end = 0;
    }
    seqSumData(vector<vector<int>>& f, double st, double en, vector<int> nam,
               string method) {
       
        summary = f;
        start = st;
        end = en;
        counts = nam;
    }
};
//**************************************************************************
DataFrame Summary::summarizeFasta(Rcpp::DataFrame& summary_report,
                                  Rcpp::IntegerVector& c) {
        vector<int> counts = Rcpp::as<std::vector<int> >(c);

        numUniques = c.size();
        total = numUniques;
        if (counts.size() != 0) {
            int initial_sum = 0;
            total = accumulate(counts.begin(), counts.end(), initial_sum);
        }else {
            counts.resize(numUniques, 1);
        }

        // columns by rows
        vector<vector<int>> report(summary_report.size());
        for (int i = 0; i < summary_report.size(); i++) {
            report[i] = Rcpp::as<std::vector<int> >(summary_report[i]);
        }

        createThreadsFasta(report, counts);

         // can't pass long long, convert to string
         vector<long long> numSeqs = getDefaults();
         vector<string> wrappedNumSeqs(numSeqs.size()+1, "");
         for (int i = 0; i < numSeqs.size(); i++) {
             wrappedNumSeqs[i] = toString(numSeqs[i]);
         }

         DataFrame df = DataFrame::create(Named("starts") = getStart(),
                                         _["ends"] = getEnd(),
                                         _["nbases"] = getLength(),
                                         _["ambigs"] = getAmbig(),
                                         _["polymers"] = getHomop(),
                                         _["numns"] = getNumNs(),
                                         _["numseqs"] = wrappedNumSeqs);
        return (df);
}
//**************************************************************************

void driverSummarize(seqSumData* params) {
    // calc lengths
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisSeqLength = params->summary[0][i];
        auto it = params->seqLength.find(thisSeqLength);
        if (it == params->seqLength.end()) {
            params->seqLength[thisSeqLength] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // starts
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisStartPosition = params->summary[1][i];
        auto it = params->startPosition.find(thisStartPosition);
        if (it == params->startPosition.end()) {
            params->startPosition[thisStartPosition] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // ends
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisEndPosition = params->summary[2][i];
        auto it = params->endPosition.find(thisEndPosition);
        if (it == params->endPosition.end()) {
            params->endPosition[thisEndPosition] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // ambigs
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisAmbig = params->summary[3][i];
        auto it = params->ambigBases.find(thisAmbig);
        if (it == params->ambigBases.end()) {
            params->ambigBases[thisAmbig] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // polymers
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisHomoP = params->summary[4][i];
        auto it = params->longHomoPolymer.find(thisHomoP);
        if (it == params->longHomoPolymer.end()) {
            params->longHomoPolymer[thisHomoP] = num; }
        else { it->second += num; } //add counts

    }

    RcppThread::checkUserInterrupt();

    // numNs
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisNumNs = params->summary[5][i];
        auto it = params->numNs.find(thisNumNs);
        if (it == params->numNs.end()) {
            params->numNs[thisNumNs] = num; }
        else { it->second += num; } //add counts
    }
}
//**************************************************************************
void Summary::createThreadsFasta(vector<vector<int>>& summary,
                                 vector<int>& counts) {
         //divide reads between processors
        Utils util;
        vector<pieceOfWork> startEndIndexes = util.divideWork(counts.size(),
                                                              processors);

        vector<RcppThread::Thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {

            seqSumData* dataBundle = new seqSumData(summary,
                                                    startEndIndexes[i+1].start,
                                                    startEndIndexes[i+1].end,
                                                    counts, "fasta");
            data.push_back(dataBundle);

            workerThreads.push_back(new RcppThread::Thread(driverSummarize,
                                                           dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(summary, startEndIndexes[0].start,
                                                startEndIndexes[0].end,
                                                counts, "fasta");
        driverSummarize(dataBundle);
        startPosition = dataBundle->startPosition;
        endPosition = dataBundle->endPosition;
        seqLength = dataBundle->seqLength;
        ambigBases = dataBundle->ambigBases;
        longHomoPolymer = dataBundle->longHomoPolymer;
        numNs = dataBundle->numNs;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();

            for (auto it = data[i]->startPosition.begin();
                 it != data[i]->startPosition.end(); it++) {
                auto itMain = startPosition.find(it->first);
                if (itMain == startPosition.end()) { //newValue
                    startPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->endPosition.begin();
                 it != data[i]->endPosition.end(); it++) {
                auto itMain = endPosition.find(it->first);
                if (itMain == endPosition.end()) { //newValue
                    endPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->seqLength.begin();
                 it != data[i]->seqLength.end(); it++)		{
                auto itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->ambigBases.begin();
                 it != data[i]->ambigBases.end(); it++)		{
                auto itMain = ambigBases.find(it->first);
                if (itMain == ambigBases.end()) { //newValue
                    ambigBases[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->longHomoPolymer.begin();
                 it != data[i]->longHomoPolymer.end(); it++)		{
                auto itMain = longHomoPolymer.find(it->first);
                if (itMain == longHomoPolymer.end()) { //newValue
                    longHomoPolymer[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            for (auto it = data[i]->numNs.begin();
                 it != data[i]->numNs.end(); it++)		{
                auto itMain = numNs.find(it->first);
                if (itMain == numNs.end()) { //newValue
                    numNs[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }
}
//**************************************************************************
