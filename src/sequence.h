#ifndef SEQUENCE_H
#define SEQUENCE_H

/*
 *  sequence.h
 *
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2024 Patrick D. Schloss. All rights reserved.
 *
 *  Data Structure for a fasta sequence.
 *
 *	A sequence object has three components:
 *	i) an accession number / name,
 *	ii) the unaligned primary sequence,
 *	iii) a pairwise aligned sequence, and
 *	iv) a sequence that is aligned to a reference alignment.
 *
 */

#include "mothur.h"
#include "utils.h"

/*****************************************************************************/

class Sequence {

public:

	Sequence();
	Sequence(string, string);
    ~Sequence() = default;

	void setAligned(string);
	void setComment(string c)   { comment = c;              }
	void setName(string);
	void setPairwise(string s)  { pairwise = s;             }
	void setUnaligned(string);

	string getAligned()         { return aligned;           }
	string getName()            { return name;              }
	string getPairwise()        { return pairwise;          }
    string getUnaligned()       { return unaligned;         }
    string getComment()         { return comment;           }
	int getNumBases()           { return numBases;          }
	int getAlignLength()        { return alignmentLength;   }

    bool isAligned();
	string getInlineSeq();
    int getNumNs();
	int getStartPos();
	int getEndPos();
    void reverseComplement();
    void trim(int);
	void padToPos(int);
	void padFromPos(int);
    void filterToPos(int);   //any character BEFORE the pos is changed to '.'
    void filterFromPos(int); //any character AFTER the pos is changed to '.'

	int getAmbigBases();
	void removeAmbigBases();
	int getLongHomoPolymer();
    string convert2ints();

protected:

    string name;
    string unaligned;
    string aligned;
    string pairwise;
    string comment;
    int numBases;
    int alignmentLength;
    int longHomoPolymer;
    int ambigBases;
    int startPos, endPos;

    void initialize();
};

/*****************************************************************************/

#endif
