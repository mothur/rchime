/*
 *  The AlignChime class is based on alignchime.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 2/26/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "alignchime.h"

/******************************************************************************/
// dn = 1.4, xn = 8, xa = 1
double AlignChimes::GetScore2(double Y, double N, double A) {
	double denominator = (opt->getXn()*(N + opt->getDn()) + opt->getXa()*A);
	return (Y/denominator);
}
/******************************************************************************/
ChimeHit2 AlignChimes::alignChime(const SeqData &QSD, const SeqData &ASD, const SeqData &B_SD,
 const string &PathQA, const string &PathQB) {

	string Q3;
	string A3;
	string B3;
	Make3way m3w;
	m3w.Make3Way(QSD, ASD, B_SD, PathQA, PathQB, Q3, A3, B3);

	ChimeHit2 Hit;
	Hit.QLabel = QSD.getName();
	Hit.AbQ = QSD.getAbund();
	Hit.AbA = ASD.getAbund();
	Hit.AbB = B_SD.getAbund();
	const unsigned ColCount = SIZE(Q3);

    // Discard terminal gaps
	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;
	for (unsigned Col = 2; Col + 2 < ColCount; ++Col) {
		char q = Q3[Col];
		char a = A3[Col];
		char b = B3[Col];

		if (isacgt(q) && isacgt(a) && isacgt(b)) {
			if (ColLo == UINT_MAX) { ColLo = Col; }
			ColHi = Col;
		}
	}

	if (ColLo == UINT_MAX) { return Hit; }

	unsigned QPos = 0;
	unsigned APos = 0;
	unsigned BPos = 0;

	vector<unsigned> ColToQPos(ColLo, UINT_MAX);
	vector<unsigned> AccumCount(ColLo, UINT_MAX);
	vector<unsigned> AccumSameA(ColLo, UINT_MAX);
	vector<unsigned> AccumSameB(ColLo, UINT_MAX);
	vector<unsigned> AccumForA(ColLo, UINT_MAX);
	vector<unsigned> AccumForB(ColLo, UINT_MAX);
	vector<unsigned> AccumAbstain(ColLo, UINT_MAX);
	vector<unsigned> AccumAgainst(ColLo, UINT_MAX);

	unsigned SumSameA = 0;
	unsigned SumSameB = 0;
	unsigned SumSameAB = 0;
	unsigned Sum = 0;
	unsigned SumForA = 0;
	unsigned SumForB = 0;
	unsigned SumAbstain = 0;
	unsigned SumAgainst = 0;
	
	for (unsigned Col = ColLo; Col <= ColHi; ++Col) {
		char q = Q3[Col];
		char a = A3[Col];
		char b = B3[Col];

		if (isacgt(q) && isacgt(a) && isacgt(b)) {
			if (q == a)              { ++SumSameA;   }
			if (q == b)              { ++SumSameB;   }
			if (a == b)              { ++SumSameAB;  }
			if (q == a && q != b)    { ++SumForA;    }
			if (q == b && q != a)    { ++SumForB;    }
			if (a == b && q != a)    { ++SumAgainst; }
			if (q != a && q != b)    { ++SumAbstain; }
			++Sum;
		}

		ColToQPos.push_back(QPos);
		AccumSameA.push_back(SumSameA);
		AccumSameB.push_back(SumSameB);
		AccumCount.push_back(Sum);
		AccumForA.push_back(SumForA);
		AccumForB.push_back(SumForB);
		AccumAbstain.push_back(SumAbstain);
		AccumAgainst.push_back(SumAgainst);

		if (q != '-') { ++QPos; }
		if (a != '-') { ++APos; }
		if (b != '-') { ++BPos; }
	}

	double IdQA = double(SumSameA)/Sum;
	double IdQB = double(SumSameB)/Sum;
	double IdAB = double(SumSameAB)/Sum;
	double MaxId = max(IdQA, IdQB);

	unsigned BestXLo = UINT_MAX;
	unsigned BestXHi = UINT_MAX;
	double BestDiv = 0.0;
	double BestIdQM = 0.0;
	double BestScore = 0.0;

// Find range of cols BestXLo..BestXHi that maximizes score
	bool FirstA = false;

// NOTE: Must be < ColHi not <= because use Col+1 below
	for (unsigned Col = ColLo; Col < ColHi; ++Col) {
		
		unsigned SameAL = AccumSameA[Col];
		unsigned SameBL = AccumSameB[Col];
		unsigned SameAR = SumSameA - AccumSameA[Col];
		unsigned SameBR = SumSameB - AccumSameB[Col];

		double IdAB = double(SameAL + SameBR)/Sum;
		double IdBA = double(SameBL + SameAR)/Sum;

		unsigned ForAL = AccumForA[Col];
		unsigned ForBL = AccumForB[Col];
		unsigned ForAR = SumForA - AccumForA[Col+1];
		unsigned ForBR = SumForB - AccumForB[Col+1];
		unsigned AbL = AccumAbstain[Col];
		unsigned AbR = SumAbstain - AccumAbstain[Col+1];

		double ScoreAB = GetScore2(ForAL, ForBL, AbL)*GetScore2(ForBR, ForAR, AbR);
		double ScoreBA = GetScore2(ForBL, ForAL, AbL)*GetScore2(ForAR, ForBR, AbR);

		double DivAB = IdAB/MaxId;
		double DivBA = IdBA/MaxId;
		double MaxDiv = max(DivAB, DivBA);
		double MaxScore = max(ScoreAB, ScoreBA);

		if (MaxScore > BestScore) {
			BestScore = MaxScore;
			BestXLo = Col;
			BestXHi = Col;
			FirstA = (ScoreAB > ScoreBA);

			if (FirstA) { BestIdQM = IdAB; }
			else        { BestIdQM = IdBA; }

			if (MaxDiv > BestDiv) { BestDiv = MaxDiv; }

		} else if (MaxScore == BestScore) {
			BestXHi = Col;
			if (MaxDiv > BestDiv) { BestDiv = MaxDiv; }
		}
	}

	if (BestXLo == UINT_MAX) { return Hit; }

	// Find maximum region of identity within BestXLo..BestXHi
	unsigned ColXLo = (BestXLo + BestXHi)/2;
	unsigned ColXHi = ColXLo;
	unsigned SegLo = UINT_MAX;
	unsigned SegHi = UINT_MAX;

	for (unsigned Col = BestXLo; Col <= BestXHi; ++Col) {
		char q = Q3[Col];
		char a = A3[Col];
		char b = B3[Col];

		if (q == a && q == b) {
			if (SegLo == UINT_MAX) { SegLo = Col; }
			SegHi = Col;
		} else {
			unsigned SegLength = SegHi - SegLo + 1;
			unsigned BestSegLength = ColXHi - ColXLo + 1;
			if (SegLength > BestSegLength) {
				ColXLo = SegLo;
				ColXHi = SegHi;
			}
			SegLo = UINT_MAX;
			SegHi = UINT_MAX;
		}
	}
	unsigned SegLength = SegHi - SegLo + 1;
	unsigned BestSegLength = ColXHi - ColXLo + 1;
	if (SegLength > BestSegLength) {
		ColXLo = SegLo;
		ColXHi = SegHi;
	}

	QPos = 0;
	for (unsigned x = 0; x < ColCount; ++x) {
		if (x == ColXLo)       { Hit.QXLo = QPos; }
		else if (x == ColXHi)  { Hit.QXHi = QPos; break; }
		
		char q = Q3[x];
		if (q != '-') { ++QPos; }
	}

	Hit.ColXLo = ColXLo;
	Hit.ColXHi = ColXHi;
	Hit.PctIdAB = IdAB*100.0;
	Hit.PctIdQM = BestIdQM*100.0;
	Hit.Div = (BestDiv - 1.0)*100.0;

	Hit.Q3 = Q3;
	Hit.QLabel = QSD.getName();
	if (FirstA) {
		Hit.A3 = A3;
		Hit.B3 = B3;
		Hit.ALabel = ASD.getName();
		Hit.BLabel = B_SD.getName();
		Hit.PctIdQA = IdQA*100.0;
		Hit.PctIdQB = IdQB*100.0;
	} else {
		Hit.A3 = B3;
		Hit.B3 = A3;
		Hit.ALabel = B_SD.getName();
		Hit.BLabel = ASD.getName();
		Hit.PctIdQA = IdQB*100.0;
		Hit.PctIdQB = IdQA*100.0;
	}

// CS SNPs
	Hit.CS_LY = 0;
	Hit.CS_LN = 0;
	Hit.CS_RY = 0;
	Hit.CS_RN = 0;
	Hit.CS_LA = 0;
	Hit.CS_RA = 0;

	for (unsigned Col = ColLo; Col <= ColHi; ++Col) {

		char q = Q3[Col];
		char a = A3[Col];
		char b = B3[Col];
		if (q == a && q == b && a == b) {
			continue;
		}

		unsigned ngaps = 0;
		if (isgap(q)) { ++ngaps; }
		if (isgap(a)) { ++ngaps; }
		if (isgap(b)) { ++ngaps; }

		if (opt->getSkipgaps()) { if (ngaps == 3) { continue; } }
		else                    { if (ngaps == 2) { continue; } }
		
		if (!FirstA) {
			swap(a, b);
		}

		if (opt->getSkipgaps2()) {
			if (Col > 0 && (isgap(Q3[Col-1]) || isgap(A3[Col-1]) || isgap(B3[Col-1]))) {
				continue; }
			if (Col + 1 < ColCount && (isgap(Q3[Col+1]) || isgap(A3[Col+1]) || isgap(B3[Col+1]))) {
				continue; }
		}

		if (Col < ColXLo) {
			if (q == a && q != b)        { ++Hit.CS_LY; }
			else if (q == b && q != a)   { ++Hit.CS_LN; }
			else                         { ++Hit.CS_LA; }
		}else if (Col > ColXHi) {
			if (q == b && q != a)        { ++Hit.CS_RY; }
			else if (q == a && q != b)   { ++Hit.CS_RN; }
			else                         { ++Hit.CS_RA; }
		}
	}

	double ScoreL = GetScore2(Hit.CS_LY, Hit.CS_LN, Hit.CS_LA);
	double ScoreR = GetScore2(Hit.CS_RY, Hit.CS_RN, Hit.CS_RA);
	Hit.Score = ScoreL*ScoreR;

	return Hit;
}
/******************************************************************************/
