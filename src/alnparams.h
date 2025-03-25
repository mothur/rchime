#ifndef alnparams_h
#define alnparams_h

/*
 *  The AlnParams class is based on alnparams.h from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 3/6/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will XXXXX
 */

#include "uchime.h"
#include "myopts.h"
#include "mxmatrix.h"
#include "myutils.h"
#include "hsp.h"

const float OBVIOUSLY_WRONG_PENALTY = 1000.0;
static const char Alphabet[] = "ACGTU";

/******************************************************************************/
class AlnParams {

public:

	AlnParams();
	~AlnParams(){}

	//const float * const *SubstMx;
	MxFloatMatrix* SubstMx; 

	bool Nucleo;
	bool NucleoSet;

// Global internal gaps
	float OpenA;
	float OpenB;

	float ExtA;
	float ExtB;

// Global terminal gaps
	float LOpenA;
	float LOpenB;
	float ROpenA;
	float ROpenB;

	float LExtA;
	float LExtB;
	float RExtA;
	float RExtB;

private:

	// Local gaps
	float LocalOpen;
	float LocalExt;

	Options* opt;
	Utilities* util;

	//Mx<float> g_SubstMxf;
	//float **g_SubstMx;
	//MxMatrix* g_SubstMx;

    void Clear();
	void setLocal(float Open, float Ext);
	void SetMxFromCmdLine(bool Nucleo);
	void setPenalties(const string &OpenStr, const string &ExtStr);
	void InitFromCmdLine(bool Nucleo);
	void setNucSubstMx(double Match, double Mismatch);
	static void ParseGapStr(const string &s, float &QI, float &QL, float &QR,
  						float &TI, float &TL, float &TR);

};
/******************************************************************************/

#endif	// alnparams_h
