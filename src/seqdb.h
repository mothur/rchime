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
#include "myutils.h"
#include "seq.h"

/******************************************************************************/
class SeqDB {

public:

	SeqDB();
	SeqDB(vector<string>& seqNames, vector<string>& seqs,
		  vector<int>& abunds, bool descending = TRUE);
	~SeqDB();

	void addSeq(string, string, int a = 1);

	SeqData getSeqData(unsigned Id) const;

	// sequence string //GetSeq
	string getSeq(unsigned SeqIndex) const;

	// sequence name //GetLabel
	string getName(unsigned SeqIndex) const;

	// sequence abundance 
	int getAbundance(unsigned SeqIndex) const;

	// sequence length //GetSeqLength
	unsigned getSeqLength(unsigned SeqIndex) const;

	// number of sequences in SeqDB GetSeqCount
	unsigned getSeqCount() const; 
	bool isNucleo() const;

	void setIsNucleo(bool nucleo) {
		isNucleotide = nucleo;
	}
	void setIsNucleoSet(bool nucleo) {
		isNucleoSet = nucleo;
	}

	private:

	Utilities* util;
    vector<int>* abunds;
	vector<string>* names;
	vector<string>* seqs;

	unsigned numSeqs;
	bool isAligned; //m_isAligned
	bool isNucleotide; //m_IsNucleo
	bool isNucleoSet; //m_IsNucleoSet

	void SetIsNucleo();
	void sortDescending();
};
/******************************************************************************/

#endif
