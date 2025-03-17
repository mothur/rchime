/*
 *  chimera_uchime.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "chimera_uchime.h"

/******************************************************************************/
ChimeraUchime::ChimeraUchime(bool derep, int proc, bool si, bool hg, Rcpp::List options) :
    Chimera(derep, proc, si, hg, options) {
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset) {
    Rcpp::List results;

    // if hasGroups then createThreads

    // else
       // get names 

       // get sequences

       // get counts

       // run as one sample

    return results;
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset,
                                         Rcpp::Environment& reference) {
    Rcpp::List results;

    return results;
}
/******************************************************************************/
