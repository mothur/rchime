#ifndef UCHIME_OUTPUT_CONVERTER
#define UCHIME_OUTPUT_CONVERTER

/*
 *  uchime_output_converter.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

/*
   This class is will replace uchime's output to the "uchimeout" file. Rather
   than writing the chimeric sequence results to a file. This function will
   store the data and create a Rcpp::Dataframe object to return to R.

   Original output was included in writechhit.cpp written for uchime.

   uchimeout headers -> DataFrame columns
   "Score Query A B IdQM IdQA IdQB IdAB IdQT LY LN LA RY RN RA Div Y"

   column descriptions can be found here,

    1	Score	Higher score means more strongly chimeric alignment.
    2	Q	Query Name.
    3	A	Parent A Name.
    4	B	Parent B Name.
    6	IdQM	Percent identity of query and the model (M) constructed as a
                segment of A and a segment of B.
    7	IdQA	Percent identity of Q and A.
    8	IdQB	Percent identity of Q and B.
    9	IdAB	Percent identity of A and B
    10	IdQT	Percent identity of Q and T.
    11	LY	Yes votes in left segment.
    12	LN	No votes in left segment.
    13	LA	Abstain votes in left segment.
    14	RY	Yes votes in right segment.
    15	RN	No votes in right segment.
    16	RA	Abstain votes in right segment.
    17	Div	Divergence, defined as (IdQM - IdQT).
    18	YN	Y, N or ?, indicating whether the query was classified as
            chimeric (Y), not chimeric (N) or borderline case (?).
            The query is classified as chimeric if h >= threshold specified by
            -minh, Div > minimum divergence specified by ‑mindiv and the number
            of diffs ( (Y+N+A) in each segment (L and R) is greater than the
            minimum specified by -mindiffs. A query is unclassified if the
            maxh > h > minh, i.e. maxh is the maximum score for a non-chimera,
            and minh is the minimum score for a chimera; in between is
            unclassified.
*/

#include "mothur.h"
#include "mothur-r.h"


/******************************************************************************/
class UchimeOutputConverter {

public:

    UchimeOutputConverter() = default;
    ~UchimeOutputConverter() = default;

    // Score, queryName, parentA, parentB, chimericStatus
    // IdQM, IdQA, IdQB, IdAB, IdQT, Div,
    // LYs, LNs, LAs, RNs, RYs, RAs
    void addOutput(double, string, string, string, string,
                   double, double, double, double, double, double,
                   int, int, int, int, int, int);

    Rcpp::DataFrame getUchimeOutput();

private:

    // Q, A, B, Y
    vector<string> queries, Aparents, Bparents, chimericStatus;

    // Score(Higher score indicates chimeric status), Divs(IdQM - IdQT)
    // IdQM, IdQA, IdQB, IdQT, IdAB - similarity scores
    vector<double> scores, IdQMs, IdQAs, IdQBs, IdABs, IdQTs, Divs;

    // LY, LN, LA, RN, RY, RA - segment votes
    vector<int> LYs, LNs, LAs, RNs, RYs, RAs;

};
/******************************************************************************/
#endif
