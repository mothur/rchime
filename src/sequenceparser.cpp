//
//  sequenceparser.cpp
//
//
//  Created by Sarah Westcott on 4/23/25.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sequenceparser.h"


/************************************************************/
SequenceParser::SequenceParser(Rcpp::Environment& dataset, bool bysample) {

    // get names
    Rcpp::Function getNames = dataset["get_names"];
    names = Rcpp::as<vector<string> >(getNames());

    // get seqs
    Rcpp::Function getSeqs = dataset["get_seqs"];
    seqs = Rcpp::as<vector<string> >(getSeqs());

    // get abunds (by sample)
    Rcpp::Function getAbunds = dataset["get_seqs_abunds"];

    if (bysample) {
        // get groups
        Rcpp::Function getGroups = dataset["get_groups"];
        groups = Rcpp::as<vector<string> >(getGroups());

        // assign group indexes
        for (int i = 0; i < groups.size(); i++) {
            groupToIndex[groups[i]] = i;
        }

        abunds = Rcpp::as<vector<vector<int>> >(getAbunds(true));
    }else {
        vector<int> temp = Rcpp::as<vector<int> >(getAbunds());
        abunds.push_back(temp);
    }

}
/************************************************************/
// get names of sequence in groups
vector<vector<string> > SequenceParser::getNames(vector<string> g) {

    if (g.size() == 0) {
        vector<vector<string> > results;
        results.push_back(names);
        return results;
    }

    // vector of names for each group passed in
    vector<vector<string> > results(g.size());

    // loop through all names
    for (int i = 0; i < names.size(); i++) {

        // for each sample
        for (int j = 0; j < g.size(); j++) {

            int groupIndex = groupToIndex[g[j]];

            // if name has sequences present in the sample, add it
            if (abunds[i][groupIndex] != 0) {
                results[j].push_back(names[i]);
            }
        }
    }

    return results;
}
/************************************************************/
// get sequences in groups
vector<vector<string> > SequenceParser::getSeqs(vector<string> g) {

    if (g.size() == 0) {
        vector<vector<string> > results;
        results.push_back(seqs);
        return results;
    }

    // vector of seqs for each group passed in
    vector<vector<string> > results(g.size());

    // loop through all seqs
    for (int i = 0; i < seqs.size(); i++) {

        // for each sample
        for (int j = 0; j < g.size(); j++) {

            int groupIndex = groupToIndex[g[j]];

            // if name has sequences present in the sample, add it
            if (abunds[i][groupIndex] != 0) {
                results[j].push_back(seqs[i]);
            }
        }
    }

    return results;
}
/************************************************************/
// get abundances of sequence in groups
vector<vector<int> > SequenceParser::getAbunds(vector<string> g) {

    if (g.size() == 0) {
        return abunds;
    }

    // vector of names for each group passed in
    vector<vector<int> > results(g.size());

    // loop through all names
    for (int i = 0; i < abunds.size(); i++) {

        // for each sample
        for (int j = 0; j < g.size(); j++) {

            int groupIndex = groupToIndex[g[j]];

            // if name has sequences present in the sample, add it
            if (abunds[i][groupIndex] != 0) {
                results[j].push_back(abunds[i][groupIndex]);
            }
        }
    }

    return results;
}
/************************************************************/
