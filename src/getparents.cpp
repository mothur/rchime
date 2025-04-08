/*
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "getparents.h"

// default - minchuck 64 chunk 4
/******************************************************************************/
vector<unsigned> GetParents::getChunkInfo(unsigned L, unsigned &Length) {
	vector<unsigned> Los;

	if (L <= opt->getMinchunk()) {
		Length = L;
		Los.push_back(0);
		return Los;
	}

	Length = (L - 1)/opt->getChunks() + 1;
	if (Length < opt->getMinchunk()) {
		Length = opt->getMinchunk();
	}

	unsigned Lo = 0;
	for (;;) {
		if (Lo + Length >= L) {
			Lo = L - Length - 1;
			Los.push_back(Lo);
			return Los;
		}
		Los.push_back(Lo);
		Lo += Length;
	}
}
/******************************************************************************/
vector<unsigned> GetParents::getCandidateParents(SeqDB* data, const SeqData & query) {
	
	vector<unsigned> parents;
	set<unsigned> ParentIndexes;

	unsigned queryLength = query.getSeqLength();
	unsigned ChunkLength;
	vector<unsigned> ChunkLos = getChunkInfo(queryLength, ChunkLength);

	// after breaking the query seqs into chunks look for parents for each chunk
	unsigned ChunkCount = SIZE(ChunkLos);
	for (unsigned ChunkIndex = 0; ChunkIndex < ChunkCount; ++ChunkIndex) {
		unsigned Lo = ChunkLos[ChunkIndex];

		string chunk = query.getSeq().substr(Lo, ChunkLength);
		Lo += ChunkLength;

		addParents(data, chunk, ParentIndexes);
	}

	// add all potential parents that meet abundance requirements
	float cutoff = opt->getAbskew()*query.getAbund();
	for (auto p = ParentIndexes.begin(); p != ParentIndexes.end(); ++p) {
		unsigned ParentIndex = *p;

		string targetName = data->getName(ParentIndex);
		float targetAbund = data->getAbundance(ParentIndex);
		if (!(targetAbund < cutoff)) {
			parents.push_back(ParentIndex);
		}
	}

	return parents;
}
/******************************************************************************/
// uses kmers to find closest parents, adds them to parent pool
void GetParents::addParents(SeqDB* DB, const string &Query,
									 set<unsigned> &ParentIndexes) {

	const unsigned SeqCount = DB->getSeqCount();
	if (SeqCount == 0)
		return;

	vector<float> WordCounts;
	vector<unsigned> Order = USort(Query, DB, WordCounts);

	unsigned TopSeqIndex = Order[0];
	float TopWordCount = WordCounts[0];
	for (unsigned i = 0; i < SeqCount; ++i) {
		unsigned SeqIndex = Order[i];
		float WordCount = WordCounts[i];
		if (TopWordCount - WordCount > 1) { return; }
		ParentIndexes.insert(SeqIndex);
	}
}
/******************************************************************************/
vector<unsigned> GetParents::USort(const string &Query, const SeqDB* DB,
											 vector<float> &WordCounts) {
	vector<unsigned> Order;
	WordCounts.clear();

	setQueryWords(Query);

	const unsigned SeqCount = DB->getSeqCount();
	vector<orderFloatAbundance> sortedVector(SeqCount);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex) {
		SeqData potentialParent = DB->getSeqData(SeqIndex);
		float WordCount = (float) getNumWordsInCommon(potentialParent);
		sortedVector[SeqIndex].index = SeqIndex;
        sortedVector[SeqIndex].abund = WordCount;
	}

    sort(sortedVector.begin(), sortedVector.end(), compareFloatAbundance);

    vector<unsigned> order(sortedVector.size(), 0);
	WordCounts.resize(sortedVector.size(), 0);
    for (int i = 0; i < WordCounts.size(); i++) {
        order[i] = sortedVector[i].index;
        WordCounts[i] = sortedVector[i].abund;
    }
	
	return order;
}
/******************************************************************************/
unsigned GetParents::getWord(string seq, unsigned offset) {
	unsigned Word = 0;

	for (unsigned i = offset; i < offset+8; ++i) {
		unsigned Letter = g_CharToLetterNucleo[(int)seq[i]];
		Word = (Word*4) + Letter;
	}
	return Word;
}
/******************************************************************************/
void GetParents::setQueryWords(const string &Query) {

	queryHasWord.clear();
	queryHasWord.resize(wordCount, 0);

	if (Query.length() <= 8) { return; }

	const unsigned L = Query.length() - 8 + 1;
	for (unsigned i = 0; i < L; ++i) {
		unsigned Word = getWord(Query, i);
		queryHasWord[Word] = 1;
	}
}
/******************************************************************************/
unsigned GetParents::getNumWordsInCommon(const SeqData &Target) {

	if (Target.getSeqLength() <= 8) {
		return 0;
	}

	unsigned Count = 0;
	const unsigned L = Target.getSeqLength() - 8 + 1;
	string seq = Target.getSeq();
	for (unsigned i = 0; i < L; ++i) {
		unsigned Word = getWord(seq, i);
		
		if (queryHasWord[Word]) { ++Count; }
	}
	return Count;
}
/******************************************************************************/
