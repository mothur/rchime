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
#include "mothur.h"
#include "utils.h"
#include "uchime_main.h"

/******************************************************************************/
struct uchimeAbunds {
    vector<int> abundances;
    bool chimeric;

    uchimeAbunds() { abundances = nullIntVector; chimeric = true; }
    uchimeAbunds(vector<int> abunds, bool c) { abundances = abunds; chimeric = c; }
    ~uchimeAbunds() {}

};
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
    Rcpp::List createUchimeResults(vector<ChimeHit2>);

    map<string, uchimeAbunds > combineResults(chimeraData*& dataBundle,
                                    vector<chimeraData*>& data,
                                    vector<RcppThread::Thread*>& workerThreads,
                                    vector<ChimeHit2>& results);
    double Pct(double x, double y);

};
/******************************************************************************/
#endif
