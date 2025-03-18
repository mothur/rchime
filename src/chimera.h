#ifndef CHIMERA_H
#define CHIMERA_H

/*
 *  chimera.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

/* This class is a parent to ChimeraUchime, ChimeraVsearch */

#include "mothur.h"
#include "mothur-r.h"
#include "rcpp_utils.h"
#include "myopts.h"

/******************************************************************************/
struct chimeraData {
    
    // inputs
    vector<string> groups;
    vector<vector<string> > names;
    vector<vector<string> > seqs;
    vector<string> refNames;
    vector<string> refSeqs;
    vector<vector<int> > abunds;

    // outputs 
    // list of dataframes storing uchimeout for each sample
    Rcpp::List uchimeOutResults;
    // list of dataframes storing uchimealns for each sample (optional)
    Rcpp::List uchimeAlnsResults;

    chimeraData(vector<string> g, vector<vector<string> > n,
                vector<vector<string> > s, vector<vector<int> > a) {
        groups = g;
        names = n;
        seqs = s;
        abunds = a;
    }
    chimeraData(vector<string> g, vector<vector<string> > n,
                vector<vector<string> > s, vector<vector<int> > a,
                vector<string>  rn, vector<string>  rs) {
        groups = g;
        names = n;
        seqs = s;
        abunds = a;
        refNames = rn;
        refSeqs = rs;
    }
    ~chimeraData() {}
};
/******************************************************************************/
class Chimera {

public:

    // dereplicate, processors, silent
    Chimera(bool derep, int proc, bool si, bool hg, Rcpp::List options);

    virtual ~Chimera(){};

    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset) = 0;
    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset,
                                      Rcpp::Environment& reference) = 0;

protected:

    int processors;
    bool silent, dereplicate, hasGroupData;
    vector<string> groups;

    Options* opts;

    void fillData(Rcpp::Environment& dataset, string group, 
                       vector<vector<string>>& names, 
                       vector<vector<string>>& seqs,
                       vector<vector<int>>& counts);

    map<string, vector<int> > combineResults(chimeraData*& dataBundle,
                                         vector<chimeraData*>& data,
                                         vector<RcppThread::Thread*>& workerThreads,
                                         Rcpp::List& results);

private:


};
/******************************************************************************/
#endif
