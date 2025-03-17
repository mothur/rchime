/*
 *  The PotentialParents class is based on addtargets2.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 *  This class receives a piece of the query and the database and
 *   returns potential chimeric parent indexes.
 *
 */

#include "seqdb.h"
#include "uchime.h"
#include "alpha.h"

const float MAX_WORD_COUNT_DROP = 1;

/******************************************************************************/
class PotentialParents {
    public:

    PotentialParents() = default;
    ~PotentialParents() = default;

    // adds potential parents to TargetIndexes
    // this function is called repeatedly from getparents with different chucks
    // of the sequence to find a pool of close potential parents
    void addTargets(SeqDB* DB, const SeqData &Query, set<unsigned> &TargetIndexes);

    private:

    vector<unsigned> queryHasWord;

    // these functions find the common words between the query and the database
    // to find the closest potential parents
    void setQueryWords(const SeqData &Query);
    unsigned getWord(string Seq, unsigned offset);
    unsigned getUniqueWordsInCommon(const SeqData &Target);
    vector<unsigned> USort(const SeqData &Query, const SeqDB* DB,
                             vector<float> &WordCounts);
    
};
/******************************************************************************/
