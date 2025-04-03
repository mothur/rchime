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
        revComp = false;
        abund = 0;
    }
    SeqData(string seqname, string sequence, unsigned abundance, unsigned ind,
            bool rev = false) { 
        name = seqname;
        seq = sequence;
        index = ind;
        revComp = rev;
        abund = abundance;
    }
    ~SeqData() = default;

    void setSeq(string sequence)  { seq = sequence; }

    //Label
    string getName() const         { return name;         }
    string getSeq() const          { return seq;          }
    unsigned getIndex() const      { return index;        }
    unsigned getSeqLength() const  { return seq.length(); }
    int  getAbund() const          { return abund;        }
    
private:

    string name;
    string seq;
    unsigned index;
    int abund;

    // RevComp means that SeqData.Seq is reverse-complemented relative
    // to the sequence in the input file (query or db). Coordinates in
    // a hit (e.g., AlnData) will be relative to SeqData.Seq, so both
    // the sequence and the coordinates should be r.c.'d for output.
    bool revComp;
};
/******************************************************************************/

#endif // seq_h
