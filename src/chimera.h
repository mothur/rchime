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
#include "options.h"
#include "chime.h"

/******************************************************************************/
struct chimeraData {
    
    // inputs
    vector<string> groups;
    vector< vector<string> >  names;
    vector< vector<string> >  seqs;
    vector<string> refNames;
    vector<string> refSeqs;
    vector< vector<int> > abunds;
    vector<set<string> > chimeras;

    // options
    Options options;

    // outputs - one for each sample
    vector<vector<ChimeHit2> > uchimeResults;
    bool silent;
    double numSeqs;

    // denovo
    chimeraData(vector< vector<string> > n, vector< vector<string> > sq,
                vector< vector<int> > a, vector<string> g, Options o,
         bool s = true) {
        groups = g;
        silent = s;
        options = o;
        names = n;
        seqs = sq;
        abunds = a;
        if (groups.size() != 0) {
            chimeras.resize(groups.size());
        }else {
            chimeras.resize(1);
        }

        numSeqs = 0;
        for (int i = 0; i < names.size(); i++) {
            numSeqs += names[i].size();
        }
    }

    // reference
    chimeraData(vector<string> n, vector<string> sq,
        vector<int> a, vector<string> rn, vector<string> rs, Options o,
        bool s = true) {
        chimeras.resize(1);
        options = o;
        silent = s;
        names.push_back(n);
        seqs.push_back(sq);
        abunds.push_back(a);
        refSeqs = rs;
        refNames = rn;

        numSeqs = names.size();
    }
    ~chimeraData() {}
};
/******************************************************************************/
class Chimera {

public:

    // dereplicate, processors, silent, hasGroups, uchime options
    Chimera(bool derep, int proc, bool si, bool hg, Rcpp::List options);

    virtual ~Chimera(){};

    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset) = 0;
    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset,
                                      Rcpp::CharacterVector& ref_names,
                                      Rcpp::CharacterVector& ref_seqs) = 0;

protected:

    int processors;
    bool silent, dereplicate, hasGroupData;
    vector<string> groups;

    Options opts;
    void removeChimerasFromDataset(Rcpp::Environment& dataset,
                                   vector<string>& seqsToRemove, string reason);
    int getNumBases(string);
    bool isGap(char c);

private:


};
/******************************************************************************/
#endif
