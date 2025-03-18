
/*
 *  chimera.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "chimera.h"

/******************************************************************************/
Chimera::Chimera(bool derep, int proc, bool si, bool hg, Rcpp::List options) {
    dereplicate = derep;
    processors = proc;
    silent = si;
    hasGroupData = hg;

    opts = Options::getInstance();

    Rcpp::CharacterVector optionsNames = options.names();
    Rcpp_Utils utils;

     bool isAligned = true;
     if (utils.contains("abskew", optionsNames))   {
         opts->setAbskew(Rcpp::as<float>(options["abskew"]));   }
     if (utils.contains("chimealns", optionsNames))   {
         opts->setChimealns(Rcpp::as<bool>(options["chimealns"]));   }
     if (utils.contains("minh", optionsNames))   {
         opts->setMinh(Rcpp::as<float>(options["minh"]));   }
     if (utils.contains("mindiv", optionsNames))   {
         opts->setMindiv(Rcpp::as<float>(options["mindiv"]));   }
     if (utils.contains("xn", optionsNames))   {
         opts->setXn(Rcpp::as<float>(options["xn"]));   }
     if (utils.contains("dn", optionsNames))   {
         opts->setDn(Rcpp::as<float>(options["dn"]));   }
     if (utils.contains("xa", optionsNames))   {
         opts->setXa(Rcpp::as<float>(options["xa"]));   }
     if (utils.contains("chunks", optionsNames))   {
         opts->setChunks(Rcpp::as<int>(options["chunks"]));   }
     if (utils.contains("minchunk", optionsNames))   {
         opts->setMinchunk(Rcpp::as<int>(options["minchunk"]));   }
     if (utils.contains("idsmoothwindow", optionsNames))   {
         opts->setIdsmoothwindow(Rcpp::as<int>(options["idsmoothwindow"]));   }
     
     if (utils.contains("maxp", optionsNames))   {
         opts->setMaxp(Rcpp::as<int>(options["maxp"]));   }
     if (utils.contains("skipgaps", optionsNames))   {
         opts->setSkipgaps(Rcpp::as<bool>(options["skipgaps"]));   }
     if (utils.contains("skipgaps2", optionsNames))   {
         opts->setSkipgaps2(Rcpp::as<bool>(options["skipgaps2"]));   }
     if (utils.contains("minlen", optionsNames))   {
         opts->setMinlen(Rcpp::as<int>(options["minlen"]));   }
     if (utils.contains("maxlen", optionsNames))   {
         opts->setMaxlen(Rcpp::as<int>(options["maxlen"]));   }

}
/******************************************************************************/
void Chimera::fillData(Rcpp::Environment& dataset, string group, 
                       vector<vector<string>>& names, 
                       vector<vector<string>>& seqs,
                       vector<vector<int>>& counts) {

    Rcpp::Function getSeqs = dataset["get_seqs"];
    Rcpp::Function getNames = dataset["get_names"];
    Rcpp::Function getCounts = dataset["get_seqs_abunds"];

    if (group == "") {
        seqs.push_back(Rcpp::as<vector<string>>(getSeqs()));
        names.push_back(Rcpp::as<vector<string>>(getNames()));
        counts.push_back(Rcpp::as<vector<int>>(getCounts()));
    }else{
        seqs.push_back(Rcpp::as<vector<string>>(getSeqs(group)));
        names.push_back(Rcpp::as<vector<string>>(getNames(group)));
        counts.push_back(Rcpp::as<vector<int>>(getCounts(group)));
    }
}
/******************************************************************************/
map<string, vector<int> > Chimera::combineResults(chimeraData*& dataBundle,
                                         vector<chimeraData*>& data,
                                         vector<RcppThread::Thread*>& workerThreads,
                                         Rcpp::List& results) {
    
    map<string, vector<int> > seqAbunds;

    // TODO 

    return seqAbunds;
}
/******************************************************************************/
