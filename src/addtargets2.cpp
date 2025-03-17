/*
 *  The PotentialParents class is based on addtargets2.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */


#include "addtargets2.h"

/******************************************************************************/
// uses kmers to find closest parents, adds them to parent pool
void PotentialParents::addTargets(SeqDB* DB, const SeqData &Query,
									 set<unsigned> &TargetIndexes) {

	const unsigned SeqCount = DB->getSeqCount();
	if (SeqCount == 0)
		return;

	vector<float> WordCounts;
	vector<unsigned> Order = USort(Query, DB, WordCounts);

	unsigned TopSeqIndex = Order[0];
	float TopWordCount = WordCounts[TopSeqIndex];
	for (unsigned i = 0; i < SeqCount; ++i) {
		unsigned SeqIndex = Order[i];
		float WordCount = WordCounts[SeqIndex];
		if (TopWordCount - WordCount > MAX_WORD_COUNT_DROP)
			return;
		TargetIndexes.insert(SeqIndex);
	}
}
/******************************************************************************/
vector<unsigned> PotentialParents::USort(const SeqData &Query, const SeqDB* DB,
											 vector<float> &WordCounts) {
	vector<unsigned> Order;
	WordCounts.clear();

	setQueryWords(Query);

	const unsigned SeqCount = DB->getSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex) {
		SeqData Target = DB->getSeqData(SeqIndex);
		float WordCount = (float) getUniqueWordsInCommon(Target);
		WordCounts.push_back(WordCount);
	}

	vector<orderFloatAbundance> sortedVector(WordCounts.size());

    for (int i = 0; i < WordCounts.size(); i++) {
        sortedVector[i].index = i;
        sortedVector[i].abund = WordCounts[i];
    }
    sort(sortedVector.begin(), sortedVector.end(), compareFloatAbundance);

    vector<unsigned> order(WordCounts.size(), 0);
    for (int i = 0; i < WordCounts.size(); i++) {
        order[i] = sortedVector[i].index;
        WordCounts[i] = sortedVector[i].abund;
    }
	
	return order;
}
/******************************************************************************/
unsigned PotentialParents::getWord(string seq, unsigned offset) {
	unsigned Word = 0;
	unsigned opt_w = 8;
	
	for (unsigned i = offset; i < opt_w; ++i) {
		unsigned Letter = g_CharToLetterNucleo[(int)seq[i]];
		Word = (Word*4) + Letter;
	}
	return Word;
}
/******************************************************************************/
void PotentialParents::setQueryWords(const SeqData &Query) {
	unsigned opt_w = 8;

	queryHasWord.resize(wordCount, 0);

	if (Query.getSeqLength() <= opt_w) {
		return;
	}

	const unsigned L = Query.getSeqLength() - opt_w + 1;
	string seq = Query.getSeq();
	for (unsigned i = 0; i < L; ++i) {
		unsigned Word = getWord(seq, i);
		queryHasWord[Word] = 1;
	}
}
/******************************************************************************/
unsigned PotentialParents::getUniqueWordsInCommon(const SeqData &Target) {

	unsigned opt_w = 8;
	if (Target.getSeqLength() <= opt_w) {
		return 0;
	}
	unsigned Count = 0;
	const unsigned L = Target.getSeqLength() - opt_w + 1;
	string seq = Target.getSeq();
	for (unsigned i = 0; i < L; ++i) {
		unsigned Word = getWord(seq, i);
		if (queryHasWord[Word]) {
			++Count;
		}
	}
	return Count;
}
/******************************************************************************/
