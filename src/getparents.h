#ifndef GETPARENTS_H
#define GETPARENTS_H

/*
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime.h"
#include "options.h"
#include "seqdb.h"
#include "alpha.h"

/******************************************************************************/
class GetParents {
    public:

    GetParents(Options); 
    ~GetParents() = default;

    vector<unsigned> getCandidateParents(SeqDB* data, const SeqData& query);
    
    vector<unsigned> getChunkInfo(unsigned L, unsigned &Length);
    unsigned getWord(string Seq, unsigned offset);
    
    private:

    int minChunk, numChunks;
    float abskew;

    vector<unsigned> queryHasWord;
    vector<unsigned> RankParents(const string &Query, const SeqDB* DB,
                             vector<float> &WordCounts);
    void setQueryWords(const string &Query);
    unsigned getNumWordsInCommon(const SeqData &Target);

    // this function is called repeatedly from getCandidateParents with different chucks
    // of the query sequence to find a pool of close potential parents
    void addParents(SeqDB* DB, const string &Query, set<unsigned> &ParentIndexes);
    
};
/******************************************************************************/
#endif
