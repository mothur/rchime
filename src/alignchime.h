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
#include "make3way.h"

/******************************************************************************/

class AlignChimes {

    public:

    AlignChimes() { opt = Options::getInstance(); }

    ChimeHit2 alignChime(const SeqData &QSD, const SeqData &ASD, const SeqData &B_SD,
     const string &PathQA, const string &PathQB);


    private:

    Options* opt;

    double GetScore2(double Y, double N, double A);
};
/******************************************************************************/

#endif
