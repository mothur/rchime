
/*
 *  uchime_output_converter.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime_output_converter.h"

/******************************************************************************/
void UchimeOutputConverter::addOutput(double score, string query,
                                      string parentA, string parentB,
                                      string status, double qm, double qa,
                                      double qb, double ab, double qt,
                                      double div, int ly,
                                      int ln, int la, int ry, int rn, int ra) {

    // higher score indicates chimeric status more likely
    scores.push_back(score);

    // query, parents and chimeric status
    queries.push_back(query);
    Aparents.push_back(parentA);
    Bparents.push_back(parentB);
    chimericStatus.push_back(status);

    // similarity
    IdQMs.push_back(qm);
    IdQAs.push_back(qa);
    IdQBs.push_back(qb);
    IdABs.push_back(ab);
    IdQTs.push_back(qt);

    // Divs(IdQM - IdQT)
    Divs.push_back(div);

    // segment votes
    LYs.push_back(ly);
    LNs.push_back(ln);
    LAs.push_back(la);
    RNs.push_back(rn);
    RYs.push_back(ry);
    RAs.push_back(ra);
}
/******************************************************************************/
Rcpp::DataFrame UchimeOutputConverter::getUchimeOutput() {

    // Score Query A B IdQM IdQA IdQB IdAB IdQT LY LN LA RY RN RA Div Y

    // create a data.frame with the chimera report data
    Rcpp::DataFrame uchime_results = DataFrame::create(
        Named("Score") = scores,
        _("Query") = queries,
        _("A") = Aparents,
        _("B") = Bparents,
        _("IdQM") = IdQMs,
        _("IdQA") = IdQAs,
        _("IdQB") = IdQBs,
        _("IdAB") = IdABs,
        _("IdQT") = IdQTs,
        _("LY") = LYs,
        _("LN") = LNs,
        _("LA") = LAs,
        _("RY") = RYs,
        _("RN") = RNs,
        _("RA") = RAs,
        _("Div") = Divs,
        _("Y") = chimericStatus);

    return uchime_results;
}
/******************************************************************************/
