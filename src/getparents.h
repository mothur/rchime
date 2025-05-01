#ifndef GETPARENTS_H
#define GETPARENTS_H

/* This class uses a kmer search to find potential
 * parents for the query sequence.
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime.h"
#include "options.h"
#include "seqdb.h"

/******************************************************************************/
class GetParents {
    public:

    GetParents(Options); 
    ~GetParents() = default;

    vector<unsigned> getCandidateParents(SeqDB* data, const SeqData& query);
    vector<unsigned> getChunkInfo(unsigned L, unsigned &Length);
    
    private:

    int minChunk, numChunks;
    float abskew;
       
};
/******************************************************************************/
#endif
