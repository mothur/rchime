/*
 *  The Make3way class is based on make3way.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 2/26/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */


#include "make3way.h"


void Make3way::Make3Way(const SeqData &QSD, const SeqData &ASD, const SeqData &B_SD,
  						const string &PathQA, const string &PathQB, string &Q3, 
						string &A3, string &B3) {
	Q3.clear();
	A3.clear();
	B3.clear();

	string querySeq = QSD.getSeq();
	string parentASeq = ASD.getSeq();
	string parentBSeq = B_SD.getSeq();

	unsigned LQ = QSD.getSeqLength();
	unsigned LA = ASD.getSeqLength();
	unsigned LB = B_SD.getSeqLength();

	vector<unsigned> InsertCountsA(LQ+1, 0);
	unsigned QPos = 0;
	for (unsigned i = 0; i < SIZE(PathQA); ++i) {
		char c = PathQA[i];
		if (c == 'M' || c == 'D') {
			++QPos;
		}else {
			++(InsertCountsA[QPos]);
		}
	}

	vector<unsigned> InsertCountsB(LQ+1, 0);
	QPos = 0;
	for (unsigned i = 0; i < SIZE(PathQB); ++i) {

		char c = PathQB[i];
		if (c == 'M' || c == 'D') {
			++QPos;
		}else {
			++(InsertCountsB[QPos]);
		}
	}

	vector<unsigned> InsertCounts;
	for (unsigned i = 0; i <= LQ; ++i) {
		unsigned is = max(InsertCountsA[i], InsertCountsB[i]);
		InsertCounts.push_back(is);
	}

	for (unsigned i = 0; i < LQ; ++i) {
		for (unsigned k = 0; k < InsertCounts[i]; ++k) {
			Q3.push_back('-');
		}
		Q3.push_back(toupper(querySeq[i]));
	}
	for (unsigned k = 0; k < InsertCounts[LQ]; ++k) {
		Q3.push_back('-');
	}
// A
	QPos = 0;
	unsigned APos = 0;
	unsigned is = 0;
	for (unsigned i = 0; i < SIZE(PathQA); ++i) {

		char c = PathQA[i];
		if (c == 'M' || c == 'D') {
			unsigned isq = InsertCounts[QPos];
			for (unsigned i = 0; i < InsertCounts[QPos]-is; ++i) {
				A3.push_back('-');
			}
			is = 0;
			++QPos;
		}
		if (c == 'M') {
			A3.push_back(toupper(parentASeq[APos++]));
		}else if (c == 'D') {
			A3.push_back('-');
		}else if (c == 'I')  {
			++is;
			A3.push_back(toupper(parentASeq[APos++]));
		}
	}
	for (unsigned k = 0; k < InsertCounts[LQ]-is; ++k) {
		A3.push_back('-');
	}
// B
	QPos = 0;
	unsigned BPos = 0;
	is = 0;
	for (unsigned i = 0; i < SIZE(PathQB); ++i) {

		char c = PathQB[i];
		if (c == 'M' || c == 'D') {
			for (unsigned i = 0; i < InsertCounts[QPos]-is; ++i) {
				B3.push_back('-');
			}
			is = 0;
			++QPos;
		}

		if (c == 'M') {
			B3.push_back(toupper(parentBSeq[BPos++]));
		} else if (c == 'D') {
			B3.push_back('-');
		}else if (c == 'I') {
			++is;
			B3.push_back(toupper(parentBSeq[BPos++]));
		}
	}
	
	for (unsigned k = 0; k < InsertCounts[LQ]-is; ++k) {
		B3.push_back('-');
	}
}
