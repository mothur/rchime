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
ChimeraUchime::ChimeraUchime(bool derep, int proc, bool si) :
    Chimera(derep, proc, si) {
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset) {
    Rcpp::List results;

    // uchime wants unaligned sequences

    return results;
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset,
                                         Rcpp::Environment& reference) {
    Rcpp::List results;

    return results;
}
/******************************************************************************/
