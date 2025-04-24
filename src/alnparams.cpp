
#include "alnparams.h"

/******************************************************************************/
AlnParams::AlnParams() {   
	setNucSubstMx(1.0, -2.0); // opt_match, opt_mismatch

	// Local
	LocalOpen = -10.0f;
	LocalExt = -1.0f;
	
	// Global
	OpenA = OpenB = -10.0;
	
	LOpenA = LOpenB = ROpenA = ROpenB = -0.5;
	ExtA = ExtB = -1.0;
	LExtA = LExtB = RExtA = RExtB = -0.5;
}
/******************************************************************************/
void AlnParams::setNucSubstMx(double Match, double Mismatch) {
	
	unsigned N = unsigned(strlen(Alphabet));

	SubstMx.resize(256, 256);
	
	for (unsigned i = 0; i < N; ++i) {
		for (unsigned j = 0; j < N; ++j) {
			float v = float(i == j ? Match : Mismatch);

			Byte ui = (Byte) toupper(Alphabet[i]);
			Byte uj = (Byte) toupper(Alphabet[j]);
			Byte li = (Byte) tolower(ui);
			Byte lj = (Byte) tolower(uj);
			ui = (Byte) toupper(ui);
			uj = (Byte) toupper(uj);

			SubstMx.matrix[ui][uj] = v;
			SubstMx.matrix[uj][ui] = v;

			SubstMx.matrix[ui][lj] = v;
			SubstMx.matrix[uj][li] = v;

			SubstMx.matrix[li][uj] = v;
			SubstMx.matrix[lj][ui] = v;

			SubstMx.matrix[li][lj] = v;
			SubstMx.matrix[lj][li] = v;
		}
	}

	for (unsigned j = 0; j < N; ++j) {
		float v = 0.0f;

		Byte ui = (Byte) 'N';
		Byte uj = (Byte) toupper(Alphabet[j]);
		Byte li = (Byte) 'n';
		Byte lj = (Byte) tolower(uj);
		ui = (Byte) toupper(ui);
		uj = (Byte) toupper(uj);

		SubstMx.matrix[ui][uj] = v;
		SubstMx.matrix[uj][ui] = v;

		SubstMx.matrix[ui][lj] = v;
		SubstMx.matrix[uj][li] = v;

		SubstMx.matrix[li][uj] = v;
		SubstMx.matrix[lj][ui] = v;

		SubstMx.matrix[li][lj] = v;
		SubstMx.matrix[lj][li] = v;
	}
}
/******************************************************************************/
