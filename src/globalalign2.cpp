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
    //util = Utilities::getInstance(); 
	g_TBBit = new MxByteMatrix();
	g_CacheLB = 0;
}
/******************************************************************************/
GlobalAligner::~GlobalAligner() { delete alignParams; }
/******************************************************************************/
bool GlobalAligner::runViterbiFast(const SeqData &Query, const SeqData &Target, PathData &PD) {
	viterbiFast(Query.getSeq(), Query.getSeqLength(),
				 Target.getSeq(), Target.getSeqLength(), PD);
	return true;
}
/******************************************************************************/
bool GlobalAligner::globalAlign(const SeqData &Query, const SeqData &Target, string &Path) {
	PathData PD; 
	runViterbiFast(Query, Target, PD);
	Path = string(PD.Start);
	return true;
}
/******************************************************************************/
bool GlobalAligner::globalAlign(const SeqData &Query, const SeqData &Target, PathData &PD) {

	PD.Clear();
	string Path;
	bool Found = globalAlign(Query, Target, Path);
	if (!Found) {
		return false;
	}
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

	if (alength*blength > 100*1000*1000) {
		return -1.0;
	}

	allocBit(alength, blength);

	MxFloatMatrix* Mx = alignParams->SubstMx;
	float OpenA = alignParams->LOpenA;
	float ExtA = alignParams->LExtA;

	//Byte **TB = g_TBBit;
	// convert to vector<float>
	//vector<float> Mrow = g_DPRow1;
	//vector<float> Drow = g_DPRow2;

// Use Mrow[-1], so...
	//Mrow[-1] = MINUS_INFINITY;
	for (unsigned j = 0; j < g_CacheLB+3; ++j) {
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
	}

// Main loop
	float M0 = float (0);
	for (unsigned i = 0; i < alength; ++i) {

		Byte a = seqA[i];
		vector<float> MxRow = Mx->matrix[a];
		float OpenB = alignParams->LOpenB;
		float ExtB = alignParams->LExtB;
		float I0 = MINUS_INFINITY;

		vector<Byte> TBrow = g_TBBit->matrix[i];
		for (unsigned j = 0; j < blength; ++j) {
			Byte b = seqB[j];
			Byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
			// M0 = DPM[i][j]
			// I0 = DPI[i][j]
			// Drow[j] = DPD[i][j]

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

			TBrow[j] = TraceBits;
		}

	// Special case for end of Drow[]
		{
		TBrow[blength] = 0;
		float md = M0 + alignParams->ROpenB;
		Drow[blength] += alignParams->RExtB;
		
		if (md >= Drow[blength]) {
			Drow[blength] = md;
			TBrow[blength] = TRACEBITS_MD;
		}
		}

		M0 = MINUS_INFINITY;

		OpenA = alignParams->OpenA;
		ExtA = alignParams->ExtA;
	}

// Special case for last row of DPI
	vector<Byte>TBrow = g_TBBit->matrix[alength];
	float I1 = MINUS_INFINITY;

	for (unsigned j = 1; j < blength; ++j) {

		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + alignParams->ROpenA;
		I1 += alignParams->RExtA;
		
		if (mi > I1) {
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
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
void GlobalAligner::allocBit(unsigned LA, unsigned LB) {

	//g_Mx_TBBit.Alloc("TBBit", LA+1, LB+1);
	g_TBBit->resize(LA+1, LB+1);
	//g_TBBit = g_Mx_TBBit.GetData();

	if (LB > g_CacheLB) {

		//util->myfree(g_DPBuffer1, g_CacheLB);
		//util->myfree(g_DPBuffer2, g_CacheLB);
		Mrow.clear();
		Drow.clear();

		g_CacheLB = LB + 128;

		// Allow use of [-1]
		//g_DPBuffer1 = myalloc<float>(g_CacheLB+3);
		//g_DPBuffer2 = myalloc<float>(g_CacheLB+3);
		//g_DPBuffer1 = static_cast<float*>(util->mymalloc(sizeof(float) * (g_CacheLB+3)));
		//g_DPBuffer2 = static_cast<float*>(util->mymalloc(sizeof(float) * (g_CacheLB+3)));
		Mrow.resize(g_CacheLB+3, MINUS_INFINITY);
		Drow.resize(g_CacheLB+3, MINUS_INFINITY);

		//g_DPRow1 = g_DPBuffer1 + 1;
		//g_DPRow2 = g_DPBuffer2 + 1;
	}
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
