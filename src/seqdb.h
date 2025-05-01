#ifndef seqdb_h
#define seqdb_h

/*
 *
 * It is used to store the reference sequences. In denovo mode the reference
 * is filled with more abundant sequences than the query.
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

	// set createKmerDB = false for query dataset to save time / memory
	// set createKmerDB = true for reference data
	SeqDB(bool denovo, bool createKmerDB);
	SeqDB(vector<string>& seqNames, vector<string>& seqs,
		  vector<int>& abunds, bool denovo, bool createKmerDB);

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

	// if createKmerDB == true, then you can search the database for closest parents
	void addPotentialParents(string seq, set<unsigned>& parents);

	private:

    vector<int> abunds;
	vector<string> names;
	vector<string> seqs;

	bool createKmerDB;
	vector<vector<int> > kmerDB;
	vector<vector<int> > kmerDBCounts;

	unsigned numSeqs;
	bool IsDenovo;
	int maxKmer;

	void sortDescending();
	void fillKmerDB();
	void addKmerDB(string seq, int);
};
/******************************************************************************/

#endif
