#ifndef ALINECHIME_H
#define ALINECHIME_H

/*
 *  The AlignChime class is based on alignchime.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 2/26/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime.h"
#include "seqdb.h"
#include "chime.h"
//#include "make3way.h"

/******************************************************************************/

class AlignChimes {

    public:

    AlignChimes(double, double, double, bool, bool); 
    ~AlignChimes() = default;

    ChimeHit2 alignChime(const SeqData &QSD, const SeqData &ASD, const SeqData &B_SD,
     const string &PathQA, const string &PathQB);

    private:

    double xn, dn, xa;
    bool skipgaps, skipgaps2;

    double GetScore2(double Y, double N, double A);
    void make3Way(const SeqData &QSD, const SeqData &ASD, const SeqData &B_SD,
        const string &PathQA, const string &PathQB, string &Q3, 
      string &A3, string &B3);

};
/******************************************************************************/

#endif
