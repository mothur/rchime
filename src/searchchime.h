#ifndef SEARCH_CHIME_H
#define SEARCH_CHIME_H

/*
 *  searchchime.h
 *
 *  The SearchChime class is based on searchchime.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 3/3/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 *  This class will evaluate the query sequence to determine if its chimeric.
 *  The results are stored in the ChimeHit2 class.
 *
 */


#include "seqdb.h"
#include "uchime.h"
#include "chime.h"
#include "getparents.h"
#include "myopts.h"
#include "path.h"
#include "globalalign2.h"
#include "alignchime.h"
#include "alpha.h"

/******************************************************************************/
class SearchChime {

public:

    SearchChime(); 
    ~SearchChime();  

    // database, query, results
    bool searchChime(SeqDB* U, const SeqData &QSD, ChimeHit2 &Hit);

private:

    Options* opt;
    GlobalAligner* aligner;
    AlignChimes* chimeAlign;

    double getFractIdGivenPath(string A, string B, const char *Path);
};
/******************************************************************************/
#endif
