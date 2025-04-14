#ifndef seqdb_h
#define seqdb_h

/*
 *  The SeqDB class is based on seqdb.h, seqdb.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 *  This class stores the sequences for processing.
 */

#include "uchime.h"
#include "seq.h"

/******************************************************************************/
class SeqDB {

public:

    // for denovo
	SeqDB(bool denovo);
	// for by reference 
	SeqDB(vector<string>& seqNames, vector<string>& seqs,
		  vector<int>& abunds, bool denovo);
	~SeqDB();

	void addSeq(string, string, int a = 1);

	SeqData getSeqData(unsigned Id) const;

	// sequence string //GetSeq
	string getSeq(unsigned SeqIndex) const;

	// sequence name //GetLabel
	string getName(unsigned SeqIndex) const;

	// sequence abundance 
	int getAbundance(unsigned SeqIndex) const;

	// number of sequences in SeqDB GetSeqCount
	unsigned getSeqCount() const; 

	bool isDenovo() { return IsDenovo; }

	private:

    vector<int> abunds;
	vector<string> names;
	vector<string> seqs;

	unsigned numSeqs;
	bool IsDenovo; 

	void sortDescending();
};
/******************************************************************************/

#endif
