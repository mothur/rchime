
#include "seqdb.h"
#include "alpha.h"
#include <__config>

//inline int compare_doubles1(const void* a, const void* b)
//{
   // double* arg1 = (double *) a;
  //  double* arg2 = (double *) b;
  //  if( *arg1 < *arg2 ) return -1;
  //  else if( *arg1 == *arg2 ) return 0;
  ////  else return 1;
//}
//inline bool compareAbunds();//

/******************************************************************************/
SeqDB::SeqDB() {
	util = Utilities::getInstance();

	names = new std::vector<string>;
	seqs = new std::vector<string>;
	abunds = new std::vector<int>;
	numSeqs = 0;

	isAligned = false;
	isNucleotide = false;
	isNucleoSet = false;
}
/******************************************************************************/
SeqDB::SeqDB(vector<string>& seqNames, vector<string>& sequences,
			   vector<int>& ab, bool descending) {
	
	util = Utilities::getInstance();

	names = new std::vector<string>;
	seqs = new std::vector<string>;
	abunds = new std::vector<int>;

	names = &seqNames;
	seqs = &sequences;
	abunds = &ab;
	numSeqs = seqNames.size();
	
	isAligned = false;
	isNucleotide = false;
	isNucleoSet = false;

	SetIsNucleo();

	if (descending) {
		sortDescending();
	}
}
/******************************************************************************/
SeqDB::~SeqDB() {
	delete names;
	delete seqs;
	delete abunds;
}
/******************************************************************************/
void SeqDB::addSeq(string name, string seq, int abund) {
	names->push_back(name);
	seqs->push_back(seq);
	abunds->push_back(abund);
	numSeqs++;
}
/******************************************************************************/
// sequence string //GetSeq
string SeqDB::getSeq(unsigned id) const {
	return (*seqs)[id];
}
/******************************************************************************/
// sequence name //GetLabel
string SeqDB::getName(unsigned id) const {
	return (*names)[id];
}
/******************************************************************************/
int SeqDB::getAbundance(unsigned id) const {
	return (*abunds)[id];
}
/******************************************************************************/
// sequence length //GetSeqLength
unsigned SeqDB::getSeqLength(unsigned id) const {
	return (*seqs)[id].length();
}
/******************************************************************************/
// number of sequences in SeqDB GetSeqCount
unsigned SeqDB::getSeqCount() const {
	return numSeqs;
}
/******************************************************************************/
bool SeqDB::isNucleo() const {
	return isNucleotide;
}
/******************************************************************************/
SeqData SeqDB::getSeqData(unsigned id) const {
	//string seqname, string sequence, unsigned abund, unsigned ind,
     //       bool rev = false, bool isNucleo = false
	SeqData seq(getName(id), getSeq(id), getAbundance(id), id, false, isNucleotide);
	return seq;
}
/******************************************************************************/
void SeqDB::sortDescending() {

	vector<orderAbundance> sortedVector((*abunds).size());

    for (int i = 0; i < (*abunds).size(); i++) {
        sortedVector[i].index = i;
        sortedVector[i].abund = (*abunds)[i];
    }
    sort(sortedVector.begin(), sortedVector.end(), compareAbundance);

    vector<unsigned> order((*abunds).size(), 0);
    for (int i = 0; i < (*abunds).size(); i++) {
        order[i] = sortedVector[i].index;
        (*abunds)[i] = sortedVector[i].abund;
    }
    
    applyOrder(*names, order);
    applyOrder(*seqs, order);
    applyOrder(*abunds, order);
}
/******************************************************************************/
void SeqDB::SetIsNucleo() {
	
	unsigned N = 0;
	for (unsigned i = 0; i < 100; ++i) {
		unsigned randomSeqIndex = util->getRandomIndex(numSeqs-1);

		string seq = getSeq(randomSeqIndex);
		unsigned L = seq.length();
		const unsigned Pos = util->getRandomIndex(L-1);
		Byte c = seq[Pos];

		if (g_IsNucleoChar[c])
			++N;
		}
	isNucleotide = (N > 80);
	isNucleoSet = true;
}
/******************************************************************************/
