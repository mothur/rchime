#ifndef CHIMERAUCHIME_H
#define CHIMERAUCHIME_H

/*
 *  chimera_uchime.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

/* This class is inherits from Chimera */


#include "chimera.h"
#include "utils.h"
#include "uchime_main.h"

/******************************************************************************/
class ChimeraUchime : public Chimera {

public:

    // dereplicate, processors, silent
    ChimeraUchime(bool derep, int proc, bool si, bool hg, Rcpp::List options);

    Rcpp::List removeChimeras(Rcpp::Environment& dataset);
    Rcpp::List removeChimeras(Rcpp::Environment& dataset,
                                      Rcpp::Environment& reference);

private:

    vector<ChimeHit2> createProcesses(Rcpp::Environment& dataset);
    Rcpp::DataFrame createUchimeOut(vector<ChimeHit2>);
    Rcpp::DataFrame createUchimeAlns(vector<ChimeHit2>);

    map<string, vector<int> > combineResults(chimeraData*& dataBundle,
                                    vector<chimeraData*>& data,
                                    vector<RcppThread::Thread*>& workerThreads,
                                    vector<ChimeHit2>& results);

};
/******************************************************************************/
#endif
