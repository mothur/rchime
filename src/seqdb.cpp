
#include "seqdb.h"
#include "alpha.h"
#include <__config>

/******************************************************************************/
SeqDB::SeqDB() {
	numSeqs = 0;
	isAligned = false;
	//isNucleotide = false;
	//isNucleoSet = false;
}
/******************************************************************************/
SeqDB::SeqDB(vector<string>& seqNames, vector<string>& sequences,
			   vector<int>& ab, bool descending) {

	//util = Utilities::getInstance();

	// degap sequences
	for (int i = 0; i < sequences.size(); i++) {
		string temp = sequences[i];
		sequences[i] = "";
		for(int j = 0; j < temp.length(); j++) {
			if(isalpha(temp[j]))	{	sequences[i] += temp[j];	}
		}
	}

	names = seqNames;
	seqs = sequences;
	abunds = ab;
	numSeqs = seqNames.size();

	isAligned = false;
	//isNucleotide = false;
	//isNucleoSet = false;

	//SetIsNucleo();

	if (descending) {
		sortDescending();
	}
}
/******************************************************************************/
SeqDB::~SeqDB() {}
/******************************************************************************/
void SeqDB::addSeq(string name, string seq, int abund, bool degap) {
	// degap sequences

	if (degap) {
		string temp = seq;
		seq = "";
		for(int j = 0; j < temp.length(); j++) {
			if(isalpha(temp[j]))	{	seq += temp[j];	}
		}
	}

	names.push_back(name);
	seqs.push_back(seq);
	abunds.push_back(abund);
	numSeqs++;
}
/******************************************************************************/
// sequence string //GetSeq
string SeqDB::getSeq(unsigned id) const {
	return seqs[id];
}
/******************************************************************************/
// sequence name //GetLabel
string SeqDB::getName(unsigned id) const {
	return names[id];
}
/******************************************************************************/
int SeqDB::getAbundance(unsigned id) const {
	return abunds[id];
}
/******************************************************************************/
// sequence length //GetSeqLength
unsigned SeqDB::getSeqLength(unsigned id) const {
	return seqs[id].length();
}
/******************************************************************************/
// number of sequences in SeqDB GetSeqCount
unsigned SeqDB::getSeqCount() const {
	return numSeqs;
}
/******************************************************************************/
SeqData SeqDB::getSeqData(unsigned id) const {
	SeqData seq(getName(id), getSeq(id), getAbundance(id),
	 id);
	return seq;
}
/******************************************************************************/
void SeqDB::sortDescending() {

	vector<orderAbundance> sortedVector((abunds).size());

    for (int i = 0; i < (abunds).size(); i++) {
        sortedVector[i].index = i;
        sortedVector[i].abund = (abunds)[i];
    }
    sort(sortedVector.begin(), sortedVector.end(), compareAbundance);

    vector<unsigned> order((abunds).size(), 0);
    for (int i = 0; i < (abunds).size(); i++) {
        order[i] = sortedVector[i].index;
        (abunds)[i] = sortedVector[i].abund;
    }

    applyOrder(names, order);
    applyOrder(seqs, order);
}
/******************************************************************************/
