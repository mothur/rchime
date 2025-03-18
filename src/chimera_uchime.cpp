/*
 *  chimera_uchime.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "chimera_uchime.h"
#include "mothur.h"
#include "uchime_main.h"

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
        params->uchimeOutResults.push_back(uchime->runUchime(params->names[k], params->seqs[k],
                                                 params->refNames, params->refSeqs, params->abunds[k]));
        delete uchime;
    }
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset) {
    Rcpp::List results;


    if (hasGroupData) {
        
        Rcpp::Function getGroups = dataset["get_groups"];
        groups = Rcpp::as<vector<string>>(getGroups());

        // run precluster on each sample individually
        results = createProcesses(dataset);
    
    } else {
        vector<vector<string> > names;
        vector<vector<string> > seqs;
        vector<vector<int> > counts;
        fillData(dataset, "", names, seqs, counts);

        chimeraData* dataBundle = new chimeraData(groups, names, seqs,
                                                        counts);
       // run as one sample
       processUchime(dataBundle);

       // remove sequences found to be chimeric
    }

    return results;
}
/******************************************************************************/
Rcpp::List ChimeraUchime::createProcesses(Rcpp::Environment& dataset) {
    Rcpp::List results;

    Utils util;
    vector<pieceOfWork> indexes = util.divideWork(groups.size(), processors);

    vector<RcppThread::Thread*> workerThreads;
    vector<chimeraData*> data;

    //Launch worker threads
    for (int i = 0; i < processors-1; i++) {
        vector<vector<string> > names;
        vector<vector<string> > seqs;
        vector<vector<int> > counts;

        // fill seqs with this processors groups
        vector<string> theseGroups;
        for (int j = indexes[i+1].start; j < indexes[i+1].end; j++) {
            fillData(dataset, groups[j], names, seqs, counts);
            theseGroups.push_back(groups[j]);
        }
        chimeraData* dataBundle = new chimeraData(theseGroups, names, seqs,
                                                        counts);
        data.push_back(dataBundle);
        workerThreads.push_back(new RcppThread::Thread(processUchime,
                                                       dataBundle));
    }

    vector<vector<string> > names;
    vector<vector<string> > seqs;
    vector<vector<int> > counts;

    // fill seqs with this processors groups
    vector<string> theseGroups;
    for (int j = indexes[0].start; j < indexes[0].end; j++) {
        fillData(dataset, groups[j], names, seqs, counts);
        theseGroups.push_back(groups[j]);
    }

    chimeraData* dataBundle = new chimeraData(theseGroups, names,
                                                 seqs, counts);
    processUchime(dataBundle);

    if (dereplicate) {
        map<string, vector<int> > seqAbunds = combineResults(dataBundle, data,
                                                         workerThreads, results);
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
        // merge results, if a sequence is found to be chimeric in any sample
        // then remove it from all samples
        
    }

    return results;
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset,
                                         Rcpp::Environment& reference) {
    Rcpp::List results;

    return results;
}
/******************************************************************************/
