
#include "seqdb.h"
#include "alpha.h"

/******************************************************************************/
SeqDB::SeqDB(bool denovo) {
	numSeqs = 0;
	IsDenovo = denovo;
}
/******************************************************************************/
SeqDB::SeqDB(vector<string>& seqNames, vector<string>& sequences,
			   vector<int>& ab, bool denovo) {

	IsDenovo = denovo;

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

	sortDescending();
}
/******************************************************************************/
SeqDB::~SeqDB() {}
/******************************************************************************/
void SeqDB::addSeq(string name, string seq, int abund) {
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
