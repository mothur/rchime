#ifndef seq_h
#define seq_h

/*
 *  The SeqData class is based on seq.h from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 2/25/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class stores the data for sequences
 *
 */

#include "uchime.h"

/******************************************************************************/
class SeqData {

public:

    SeqData() { 
        name = "";
        seq = "";
        index = UINT_MAX;
        abund = 0;
    }
    SeqData(string seqname, string sequence, unsigned abundance, unsigned ind) { 
        name = seqname;
        seq = sequence;
        index = ind;
        abund = abundance;
    }
    ~SeqData() = default;

    string getName() const         { return name;         }
    string getSeq() const          { return seq;          }
    void setSeq(string sequence)   { seq = sequence;      }
    unsigned getIndex() const      { return index;        }
    unsigned getSeqLength() const  { return seq.length(); }
    int  getAbund() const          { return abund;        }
    
private:

    string name;
    string seq;
    unsigned index;
    int abund;
};
/******************************************************************************/

#endif // seq_h
