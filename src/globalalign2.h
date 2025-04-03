#ifndef GLOBAL_ALIGNER
#define GLOBAL_ALIGNER

/*
 *  The GlobalAligner class is based on globalalign2.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 3/10/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will fill PathData with the alignment info
 */

#include "uchime.h"
#include "seq.h"
#include "path.h"
#include "alnparams.h"
#include "myutils.h"

const Byte TRACEBITS_DM = 0x01;
const Byte TRACEBITS_IM = 0x02;
const Byte TRACEBITS_MD = 0x04;
const Byte TRACEBITS_MI = 0x08;
const Byte TRACEBITS_SM = 0x10;
const Byte TRACEBITS_UNINIT = ~0x1f;

/******************************************************************************/
class GlobalAligner {

public:

	GlobalAligner(); 
	~GlobalAligner(); 

	bool globalAlign(const SeqData &Query, const SeqData &Target, PathData &PD);
	
private:

	AlnParams* alignParams;
    MxByteMatrix* g_TBBit;

    bool globalAlign(const SeqData &Query, const SeqData &Target, string &Path);
    float viterbiFast(string seqA, unsigned alength, string seqB, unsigned blength,
            PathData &PD);
    void traceBackBit(unsigned LA, unsigned LB, char State, PathData &PD);
};
/******************************************************************************/

#endif
