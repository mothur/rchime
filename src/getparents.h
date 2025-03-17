#ifndef GETPARENTS_H
#define GETPARENTS_H

/*
 *  The GetParents class is based on getparents.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime.h"
#include "myopts.h"
#include "seqdb.h"

/******************************************************************************/
class GetParents {
    public:

    GetParents() { opt = Options::getInstance(); }
    ~GetParents() = default;

    vector<unsigned> getChunkInfo(unsigned L, unsigned &Length);
    vector<unsigned> getCandidateParents(SeqDB* data, const SeqData& query);

    private:

    Options* opt;
};
/******************************************************************************/
#endif
