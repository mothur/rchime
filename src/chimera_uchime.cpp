/*
 *  chimera_uchime.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "chimera_uchime.h"
#include "chime.h"
#include "mothur.h"
#include "uchime_main.h"




// create a data.frame with the chimera report data
    // Rcpp::DataFrame uchime_results = DataFrame::create(
    //     Named("Score") = scores,
    //     _("Query") = queries,
    //     _("A") = Aparents,
    //     _("B") = Bparents,
    //     _("IdQM") = IdQMs,
    //     _("IdQA") = IdQAs,
    //     _("IdQB") = IdQBs,
    //     _("IdAB") = IdABs,
    //     _("IdQT") = IdQTs,
    //     _("LY") = LYs,
    //     _("LN") = LNs,
    //     _("LA") = LAs,
    //     _("RY") = RYs,
    //     _("RN") = RNs,
    //     _("RA") = RAs,
    //     _("Div") = Divs,
    //     _("Y") = chimericStatus);


    
/******************************************************************************/
ChimeraUchime::ChimeraUchime(bool derep, int proc, bool si, bool hg, Rcpp::List options) :
    Chimera(derep, proc, si, hg, options) {
}
/******************************************************************************/
void processUchime(chimeraData* params) {
    
    // for each group in dataset
    for (int k = 0; k < params->seqs.size(); k++) {
        // long long numSeqs = params->seqs[k].size();

        // call uchime source code
        UchimeMain* uchime = new UchimeMain();
        params->uchimeResults.push_back(uchime->runUchime(params->names[k], params->seqs[k],
                                                params->refNames, params->refSeqs,
                                                params->abunds[k], params->chimeras[k]));                                         
        delete uchime;
    }
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset) {
    
    vector<ChimeHit2> uchimeOutputs;
    vector<string> resultsNames; 

    if (hasGroupData) {
        
        Rcpp::Function getGroups = dataset["get_groups"];
        groups = Rcpp::as<vector<string>>(getGroups());

        // run uchime on each sample individually
        uchimeOutputs = createProcesses(dataset);
    
    } else {
        chimeraData* dataBundle = new chimeraData(dataset, groups);
        // run as one sample
        processUchime(dataBundle);

        uchimeOutputs = dataBundle->uchimeResults[0];
        vector<string> chimeras = toVector(dataBundle->chimeras[0]);

        delete dataBundle;
        removeChimerasFromDataset(dataset, chimeras, "chimera_uchime");
    }

    Rcpp::List results;
    results.push_back(createUchimeOut(uchimeOutputs));
    
    resultsNames.push_back("uchimeout");
    
    if (opts->getChimealns()) {
        results.push_back(createUchimeAlns(uchimeOutputs));
        resultsNames.push_back("uchimealns");
    }
    results.attr("names") = resultsNames;

    return results;
}
/******************************************************************************/
vector<ChimeHit2> ChimeraUchime::createProcesses(Rcpp::Environment& dataset) {
    vector<ChimeHit2> results;

    Utils util;
    vector<pieceOfWork> indexes = util.divideWork(groups.size(), processors);

    vector<RcppThread::Thread*> workerThreads;
    vector<chimeraData*> data;

    //Launch worker threads
    for (int i = 0; i < processors-1; i++) {
        // fill seqs with this processors groups
        vector<string> theseGroups;
        for (int j = indexes[i+1].start; j < indexes[i+1].end; j++) {
            theseGroups.push_back(groups[j]);
        }
        chimeraData* dataBundle = new chimeraData(dataset, theseGroups);
        data.push_back(dataBundle);
        workerThreads.push_back(new RcppThread::Thread(processUchime,
                                                       dataBundle));
    }

    // fill seqs with this processors groups
    vector<string> theseGroups;
    for (int j = indexes[0].start; j < indexes[0].end; j++) {
        theseGroups.push_back(groups[j]);
    }

    chimeraData* dataBundle = new chimeraData(dataset, theseGroups);
    processUchime(dataBundle);

    map<string, vector<int> > seqAbunds = combineResults(dataBundle, data,
                                                         workerThreads, results);

    if (dereplicate) {
        // set dataset abundances
        Rcpp::Function setAbundances = dataset["set_abundances"];
        Rcpp::List abundances;
        vector<string> names;
        for (auto it = seqAbunds.begin(); it != seqAbunds.end(); it++) {
            names.push_back(it->first);
            abundances.push_back(it->second);
        }

        setAbundances(names, abundances, Rcpp::Named("reason",
                                                 "chimera_uchime"));
    }else {
        // seqAbunds just contains seqs to be removed
        vector<string> seqsToRemove = getKeys(seqAbunds);
        removeChimerasFromDataset(dataset, seqsToRemove, "chimera_uchime");
    }

    return results;
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset,
                                         Rcpp::Environment& reference) {
    Rcpp::List results = Rcpp::List::create();
    vector<string> resultsNames;

    chimeraData* dataBundle = new chimeraData(dataset, reference);

    // run as one sample
    processUchime(dataBundle);


    vector<ChimeHit2> hits = (dataBundle->uchimeResults[0]);
    vector<string> chimeras = toVector(dataBundle->chimeras[0]);
    delete dataBundle;

    removeChimerasFromDataset(dataset, chimeras, "chimera_uchime");

    results.push_back(createUchimeOut(hits));
    resultsNames.push_back("uchimeout");
    
    if (opts->getChimealns()) {
        results.push_back(createUchimeAlns(hits));
        resultsNames.push_back("uchimealns");
    }
    results.attr("names") = resultsNames;

    return results;
}
/******************************************************************************/
// this is called when dereplicate = TRUE, meaning sequence are only removed 
// from the samples they are flagged as chimeric in
map<string, vector<int> > ChimeraUchime::combineResults(chimeraData*& dataBundle,
                                         vector<chimeraData*>& data,
                                         vector<RcppThread::Thread*>& workerThreads,
                                         vector<ChimeHit2>& results) {
    
    map<string, vector<int> > seqAbunds;
    int groupIndex = 0;
    vector<vector<ChimeHit2> > uchimeOut = dataBundle->uchimeResults;
          
    // dereplicate = TRUE, guilty until proven innocent (all samples must find it to be chimeric)
    // dereplicate = FALSE, innocent until proven guilty (any sample can find it to be chimeric)
    map<string, bool> flaggedChimericInAllSamples;
    
    // collect main threads results
    for (int i = 0; i < dataBundle->groups.size(); i++) {

        for (int j = 0; j < dataBundle->names[i].size(); j++) {

            string seqName = dataBundle->names[i][j];
            auto it = seqAbunds.find(seqName);

            // create new abunds vector and fill this groups abunds
            if (it == seqAbunds.end()) {
                vector<int> abunds(groups.size(), 0);
                flaggedChimericInAllSamples[seqName] = dereplicate;

                if (dereplicate) {
                    // if not chimeric in this sample, then add abundance
                    if (dataBundle->chimeras[i].count(seqName) == 0){
                        abunds[groupIndex] = dataBundle->abunds[i][j];
                        flaggedChimericInAllSamples[seqName] = false;
                    }
                    seqAbunds[seqName] = abunds;
                }else {
                    // if chimeric in this sample, then remove from dataset
                    if (dataBundle->chimeras[i].count(seqName) != 0){
                       seqAbunds[seqName] = nullIntVector; 
                       flaggedChimericInAllSamples[seqName] = true;
                    }
                }
            }else{
                if (dereplicate) {
                    // if not chimeric in this sample, then add abundance
                    if (dataBundle->chimeras[i].count(seqName) == 0) {
                        it->second[groupIndex] = dataBundle->abunds[i][j];
                        flaggedChimericInAllSamples[seqName] = false;
                    }
                }else {
                    // if chimeric in this sample, then remove from dataset
                    if (dataBundle->chimeras[i].count(seqName) != 0){
                       it->second = nullIntVector; 
                       flaggedChimericInAllSamples[seqName] = true;
                    }
                }
            } 
        }
        groupIndex++;
    }

    delete dataBundle;

    // collect child thread results and join
    for (int i = 0; i < processors-1; i++) {

        workerThreads[i]->join();

        for (int j = 0; j < data[i]->groups.size(); j++) {

            uchimeOut.push_back(data[i]->uchimeResults[j]);

            for (int k = 0; k < data[i]->names[j].size(); k++) {

                string seqName = data[i]->names[j][k];
                auto it = seqAbunds.find(seqName);

                // create new abunds vector and fill this groups abunds
                if (it == seqAbunds.end()) {
                    
                    vector<int> abunds(groups.size(), 0);
                    flaggedChimericInAllSamples[seqName] = dereplicate;

                    if (dereplicate) {
                        // if not chimeric in this sample, then add abundance
                        if (data[i]->chimeras[i].count(seqName) == 0){
                            abunds[groupIndex] = data[i]->abunds[i][j];
                            flaggedChimericInAllSamples[seqName] = false;
                        }
                        seqAbunds[seqName] = abunds;
                    }else {
                        // if chimeric in this sample, then remove from dataset
                        if (data[i]->chimeras[i].count(seqName) != 0){
                            seqAbunds[seqName] = nullIntVector; 
                            flaggedChimericInAllSamples[seqName] = true;
                        }
                    }
                }else{
                    if (dereplicate) {
                        // if not chimeric in this sample, then add abundance
                        if (data[i]->chimeras[i].count(seqName) == 0) {
                            it->second[groupIndex] = data[i]->abunds[i][j];
                            flaggedChimericInAllSamples[seqName] = false;
                        }
                    }else {
                        // if chimeric in this sample, then remove from dataset
                        if (data[i]->chimeras[i].count(seqName) != 0){
                            it->second = nullIntVector; 
                            flaggedChimericInAllSamples[seqName] = true;
                        }
                    }
                } 
            }
            groupIndex++;
        }

        delete data[i];
        delete workerThreads[i];
    }

    // merge uchime results using flagged chimeras
    // the result will have 1 entry in hits for each sequence in dataset

    // if all samples find the query to be chimeric, choose first chimeric result
    // else pick first non-chimeric result.

    set<string> resolved; 
    for (int i = 0; i < uchimeOut.size(); i++) {
        for (int j = 0; j < uchimeOut[i].size(); j++) {

            string seqName = uchimeOut[i][j].QLabel;
            auto it = resolved.find(seqName);

            // if we haven't recorded this sequence
            if (it == resolved.end()) {

                // sequence was found to be chimeric in all samples
                if (flaggedChimericInAllSamples[seqName]) {
                        results.push_back(uchimeOut[i][j]);
                        resolved.insert(seqName);
                }else {
                    // pick first non-chimeric results
                    if (!uchimeOut[i][j].Accept()) {
                        results.push_back(uchimeOut[i][j]);
                        resolved.insert(seqName);
                    }
                }
            }
        }
    }

    return seqAbunds;
}
/******************************************************************************/
Rcpp::DataFrame ChimeraUchime::createUchimeOut(vector<ChimeHit2>) {

    Rcpp::DataFrame output;


    return output;
    
}
/******************************************************************************/
Rcpp::DataFrame ChimeraUchime::createUchimeAlns(vector<ChimeHit2>) {
    
    Rcpp::DataFrame output;


    return output;
    
}
/******************************************************************************/
