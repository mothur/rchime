/*
 *  The GlobalAligner class is based on globalalign2.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 3/10/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will fill PathData with the alignment info
 */

#include "globalalign2.h"
#include "mxmatrix.h"

/******************************************************************************/
GlobalAligner::GlobalAligner() { 
	alignParams = new AlnParams();
	g_TBBit = new MxByteMatrix();
}
/******************************************************************************/
GlobalAligner::~GlobalAligner() { delete alignParams; }
/******************************************************************************/
bool GlobalAligner::globalAlign(const SeqData &Query, const SeqData &Target, PathData &PD) {

	PathData tempPath; 
	viterbiFast(Query.getSeq(), Query.getSeqLength(),
				 Target.getSeq(), Target.getSeqLength(), tempPath);

	string Path = string(tempPath.Start);

	unsigned n = SIZE(Path);
	PD.Alloc(n+1);
	memcpy(PD.Front, Path.c_str(), n);
	PD.Start = PD.Front;
	PD.Start[n] = 0;
	return true;
}
/******************************************************************************/
float GlobalAligner::viterbiFast(string seqA, unsigned alength,
								 string seqB, unsigned blength, PathData &PD) {

	if (alength*blength > 100*1000*1000) { return -1.0; }

	//allocBit(alength, blength);
	g_TBBit->resize(alength+1, blength+1);
	unsigned g_CacheLB = blength + 128;
	vector<float> Mrow(g_CacheLB+3, MINUS_INFINITY);
    vector<float> Drow(g_CacheLB+3, MINUS_INFINITY);

	MxFloatMatrix* Mx = alignParams->SubstMx;
	float OpenA = alignParams->LOpenA;
	float ExtA = alignParams->LExtA;

// Main loop
	float M0 = float (0);
	for (unsigned i = 0; i < alength; ++i) {

		Byte a = seqA[i];
		vector<float> MxRow = Mx->matrix[a];
		float OpenB = alignParams->LOpenB;
		float ExtB = alignParams->LExtB;
		float I0 = MINUS_INFINITY;

		for (unsigned j = 0; j < blength; ++j) {
			Byte b = seqB[j];
			Byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{

			float xM = M0;
			
			if (Drow[j] > xM) {
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
			}

			if (I0 > xM) {
				xM = I0;
				TraceBits = TRACEBITS_IM;
			}

			M0 = Mrow[j];
			Mrow[j] = xM + MxRow[b];
			}

		// DELETE
			{
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
			
			if (md >= Drow[j]) {
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
			}
			}

		// INSERT
			{
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
			
			if (mi >= I0) {
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
			}
			}

			OpenB = alignParams->OpenB;
			ExtB = alignParams->ExtB;

			g_TBBit->matrix[i][j] = TraceBits;
		}

	// Special case for end of Drow[]
		{
		g_TBBit->matrix[i][blength] = 0;
		float md = M0 + alignParams->ROpenB;
		Drow[blength] += alignParams->RExtB;
		
		if (md >= Drow[blength]) {
			Drow[blength] = md;
			g_TBBit->matrix[i][blength] = TRACEBITS_MD;
		}
		}

		M0 = MINUS_INFINITY;

		OpenA = alignParams->OpenA;
		ExtA = alignParams->ExtA;
	}

	// Special case for last row of DPI
	float I1 = MINUS_INFINITY;

	for (unsigned j = 1; j < blength; ++j) {

		g_TBBit->matrix[alength][j] = 0;
		float mi = Mrow[int(j)-1] + alignParams->ROpenA;
		I1 += alignParams->RExtA;
		
		if (mi > I1) {
			I1 = mi;
			g_TBBit->matrix[alength][j] = TRACEBITS_MI;
		}
		
	}

	float FinalM = Mrow[blength-1];
	float FinalD = Drow[blength];
	float FinalI = I1;
	float Score = FinalM;

	Byte State = 'M';
	if (FinalD > Score) {
		Score = FinalD;
		State = 'D';
	}
	if (FinalI > Score) {
		Score = FinalI;
		State = 'I';
	}

	traceBackBit(alength, blength, State, PD);

	return Score;
}
/******************************************************************************/
void GlobalAligner::traceBackBit(unsigned alength, unsigned blength, char State, PathData &PD) {

	PD.Alloc(alength+blength);

	char *PathPtr = PD.Back;
	*PathPtr = 0;

	size_t i = alength;
	size_t j = blength;
	for (;;) {

		if (i == 0 && j == 0) {
			break;
		}

		--PathPtr;
		*PathPtr = State;

		Byte t;
		switch (State) {
		case 'M':
			t = g_TBBit->matrix[i-1][j-1];
			if (t & TRACEBITS_DM) {
				State = 'D';
			} else if (t & TRACEBITS_IM) {
				State = 'I';
			}else {
				State = 'M';
			}
			--i;
			--j;
			break;
		case 'D':
			
			t = g_TBBit->matrix[i-1][j];
			if (t & TRACEBITS_MD) {
				State = 'M';
			}else{
				State = 'D';
			}
			--i;
			break;

		case 'I':
			t = g_TBBit->matrix[i][j-1];
			if (t & TRACEBITS_MI) {
				State = 'M';
			}else {
				State = 'I';
			}
			--j;
			break;

		default:
			break;
		}
	}
	PD.Start = PathPtr;
}
/******************************************************************************/
