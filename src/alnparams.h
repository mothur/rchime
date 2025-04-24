#ifndef alnparams_h
#define alnparams_h

/*
 *  The AlnParams class is based on alnparams.h from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 3/6/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class is used when aligning the query to the parents.
 */

#include "uchime.h"
#include "mxmatrix.h"

static const char Alphabet[] = "ACGTU";

/******************************************************************************/
class AlnParams {

public:

	AlnParams();
	~AlnParams(){}

	MxFloatMatrix SubstMx; 

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

	void setNucSubstMx(double Match, double Mismatch);
};
/******************************************************************************/

#endif	// alnparams_h
