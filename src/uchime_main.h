#ifndef UCHIME_MAIN_H
#define UCHIME_MAIN_H

/*
 *  uchime_main.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "myutils.h"
#include "chime.h"
#include "seqdb.h"
#include "searchchime.h"
#include "mothur-r.h"


/* This class is the access point to the uchime source */
/******************************************************************************/
class UchimeMain {

public:

    UchimeMain();   
    ~UchimeMain();  

    // returns list of dataframes with results
    vector<ChimeHit2> runUchime(vector<string> names, vector<string> seqs,
                         vector<string> refNames, vector<string> refSeqs,
                         vector<int> abunds, set<string>& chims);

private:

   Utilities* util;
   bool uchimeDeNovo;
};
/******************************************************************************/
#endif
