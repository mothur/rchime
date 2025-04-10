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

    // denovo
    chimeraData(Rcpp::Environment& d, vector<string> g, Options o,
         bool s = true) {
        groups = g;
        silent = s;
        options = o;

        if (groups.size() != 0) {
            for (int i = 0; i < groups.size(); i++) {
                fillData(d, groups[i], "dataset");
            }
            chimeras.resize(groups.size());
        }else{
            fillData(d, "", "dataset");
            chimeras.resize(1);
        }
    }

    // reference
    chimeraData(Rcpp::Environment& d, Rcpp::Environment& reference, Options o) {
        fillData(d, "", "dataset");
        fillData(reference, "", "reference");
        chimeras.resize(1);
        options = o;
    }
    ~chimeraData() {}

    void fillData(Rcpp::Environment& dataset, string group, string mode){

        Rcpp::Function getNames = dataset["get_names"];
        Rcpp::Function getSeqs = dataset["get_seqs"];

        if (mode == "dataset") {
            
            Rcpp::Function getCounts = dataset["get_seqs_abunds"];

            if (group == "") {
                seqs.push_back(Rcpp::as<vector<string> >(getSeqs()));
                names.push_back(Rcpp::as<vector<string> >(getNames()));
                abunds.push_back(Rcpp::as<vector<int> >(getCounts()));
            }else{
                seqs.push_back(Rcpp::as<vector<string> >(getSeqs(group)));
                names.push_back(Rcpp::as<vector<string> >(getNames(group)));
                abunds.push_back(Rcpp::as<vector<int> >(getCounts(group)));
            }
        }else {
            refNames = Rcpp::as<vector<string> >(getNames());
            refSeqs = Rcpp::as<vector<string> >(getSeqs());
        }
    }
    
};
/******************************************************************************/
class Chimera {

public:

    // dereplicate, processors, silent, hasGroups, uchime options
    Chimera(bool derep, int proc, bool si, bool hg, Rcpp::List options);

    virtual ~Chimera(){ delete opts; };

    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset) = 0;
    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset,
                                      Rcpp::Environment& reference) = 0;

protected:

    int processors;
    bool silent, dereplicate, hasGroupData;
    vector<string> groups;

    Options* opts;
    void removeChimerasFromDataset(Rcpp::Environment& dataset,
                                   vector<string>& seqsToRemove, string reason);
    int getNumBases(string);
    bool isGap(char c);

private:


};
/******************************************************************************/
#endif
