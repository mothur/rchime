
#include "seqdb.h"
#include "kmer.h"

/******************************************************************************/
SeqDB::SeqDB(bool denovo, bool cKmerDB) {
	numSeqs = 0;
	IsDenovo = denovo;
	createKmerDB = cKmerDB;

	if (createKmerDB) {
		int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096,
						    16384, 65536, 262144, 1048576, 4194304,
							16777216, 67108864 };

		maxKmer = power4s[8];
		kmerDB.resize(maxKmer+1);
		kmerDBCounts.resize(maxKmer+1);
	}
}
/******************************************************************************/
SeqDB::SeqDB(vector<string>& seqNames, vector<string>& sequences,
			   vector<int>& ab, bool denovo, bool cKmerDB) {

	IsDenovo = denovo;
	createKmerDB = cKmerDB;

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

	if (createKmerDB) {

		int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096,
							16384, 65536, 262144, 1048576, 4194304,
							16777216, 67108864 };

		maxKmer = power4s[8];
		kmerDB.resize(maxKmer+1);
		kmerDBCounts.resize(maxKmer+1);

		fillKmerDB();
	}
}
/******************************************************************************/
SeqDB::~SeqDB() {}
/******************************************************************************/
void SeqDB::addSeq(string name, string seq, int abund) {
	names.push_back(name);
	seqs.push_back(seq);
	abunds.push_back(abund);

	if (createKmerDB) {
		addKmerDB(seq, numSeqs);
	}
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
void SeqDB::addPotentialParents(string seq, set<unsigned>& parents) {

	if (!createKmerDB) { return; }

	vector<int> matches(numSeqs, 0);
	vector<bool> timesKmerFound(kmerDB.size()+1, false);

	Kmer kmer(8);
	int numKmers = seq.length() - 8 + 1;

	for(int i=0;i<numKmers;i++){
		int kmerNumber = kmer.getKmerNumber(seq, i);
		if(!timesKmerFound[kmerNumber]){
			for(int j=0;j<kmerDB[kmerNumber].size();j++){
				//matches[kmerDB[kmerNumber][j]]++;
				matches[kmerDB[kmerNumber][j]] += kmerDBCounts[kmerNumber][j];
			}
		}
		timesKmerFound[kmerNumber] = true;
	}

	// sort by most common kmers
	vector<orderFloatAbundance> sortedVector(numSeqs);
	for (unsigned i = 0; i < numSeqs; ++i) {
		sortedVector[i].index = i;
        sortedVector[i].abund = matches[i];
	}

    sort(sortedVector.begin(), sortedVector.end(), compareFloatAbundance);

	int TopWordCount = sortedVector[0].abund;
	for (unsigned i = 0; i < numSeqs; ++i) {
		float WordCount = sortedVector[i].abund;
		if (TopWordCount - WordCount > 1) { return; }
		parents.insert(sortedVector[i].index);
	}
}
/******************************************************************************/
void SeqDB::addKmerDB(string seq, int index) {
	Kmer kmer(8);

	int numKmers = seq.length() - 8 + 1;

	vector<int> seenBefore(maxKmer+1,0);
	for(int j=0;j<numKmers;j++){
		int kmerNumber = kmer.getKmerNumber(seq, j);
		// only add index once
		if(seenBefore[kmerNumber] == 0){
			kmerDB[kmerNumber].push_back(index);
			kmerDBCounts[kmerNumber].push_back(1);
		}else {
			// increment the number of times seen
			kmerDBCounts[kmerNumber][kmerDB[kmerNumber].size()-1]++;
		}
		seenBefore[kmerNumber] = 1;
	}
}
/******************************************************************************/
void SeqDB::fillKmerDB() {
	for (int i = 0; i < seqs.size(); i++) {
		addKmerDB(seqs[i], i);
	}
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
