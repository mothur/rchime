
#include "alnparams.h"
#include "hsp.h"

/******************************************************************************/
AlnParams::AlnParams() {   
	opt = Options::getInstance();
	util = Utilities::getInstance();
	InitFromCmdLine();
}
/******************************************************************************/
void AlnParams::Clear() {
	LocalOpen = OBVIOUSLY_WRONG_PENALTY;
	LocalExt = OBVIOUSLY_WRONG_PENALTY;
	OpenA = OBVIOUSLY_WRONG_PENALTY;
	OpenB = OBVIOUSLY_WRONG_PENALTY;
	ExtA = OBVIOUSLY_WRONG_PENALTY;
	ExtB = OBVIOUSLY_WRONG_PENALTY;
	LOpenA = OBVIOUSLY_WRONG_PENALTY;
	LOpenB = OBVIOUSLY_WRONG_PENALTY;
	ROpenA = OBVIOUSLY_WRONG_PENALTY;
	ROpenB = OBVIOUSLY_WRONG_PENALTY;
	LExtA = OBVIOUSLY_WRONG_PENALTY;
	LExtB = OBVIOUSLY_WRONG_PENALTY;
	RExtA = OBVIOUSLY_WRONG_PENALTY;
	RExtB = OBVIOUSLY_WRONG_PENALTY;
}
/******************************************************************************/
void AlnParams::setLocal(float Open, float Ext) {
	LocalOpen = Open;
	LocalExt = Ext;
}
/******************************************************************************/
/***
Open/Ext format string is one or more:
	[<flag><flag>...]<value>

Value is (positive) penalty or * (disabled).
Flag is:
	Q		Query.
	T		Target sequence.
	I		Internal gaps (defafault internal and terminal).
	E		End gaps (default internal and terminal).
	L		Left end.
	R		Right end.
***/
/******************************************************************************/
void AlnParams::ParseGapStr(const string &s, float &QI, float &QL, float &QR,
  						float &TI, float &TL, float &TR) {
	if (s.empty()) {
		return;
	}

	bool Q = false;
	bool T = false;
	bool I = false;
	bool E = false;
	bool L = false;
	bool R = false;

	const unsigned K = SIZE(s);
	unsigned Dec = 0;
	float Value = FLT_MAX;

	for (unsigned i = 0; i <= K; ++i) {
		char c = s.c_str()[i];
		if (c == 0 || c == '/') {
			if (!Q && !T && !I && !E && !L && !R) {
				Q = true;
				T = true;
				L = true;
				R = true;
				I = true;
			}

			if (!E && !I && !L && !R) {
				E = false;
				I = true;
				L = true;
				R = true;
			}

			if (E) {
				L = true;
				R = true;
			}

			if (!Q && !T) {
				Q = true;
				T = true;
			}

			if (Q && L)
				QL = -Value;
			if (Q && R)
				QR = -Value;
			if (Q && I)
				QI = -Value;
			if (T && L)
				TL = -Value;
			if (T && R)
				TR = -Value;
			if (T && I)
				TI = -Value;
			
			Value = FLT_MAX;
			Dec = 0;
			Q = false;
			T = false;
			I = false;
			E = false;
			L = false;
			R = false;
			}
		else if (c == '*') {
			Value = -MINUS_INFINITY;
		} else if (isdigit(c)) {
			if (Value == FLT_MAX) {
				Value = 0.0;
			}
			if (Dec > 0) {
				Dec *= 10;
				Value += float(c - '0')/Dec;
			} else {
				Value = Value*10 + (c - '0');
			}
		} else if (c == '.') {
			Dec = 1;
		} else {
			switch (c) {
			case 'Q':
				Q = true;
				break;
			case 'T':
				T = true;
				break;
			case 'I':
				I = true;
				break;
			case 'L':
				L = true;
				break;
			case 'R':
				R = true;
				break;
			case 'E':
				E = true;
				break;
			default:
				break;
			}
		}
	}
}
/******************************************************************************/
void AlnParams::setPenalties(const string &OpenStr, const string &ExtStr) {
	ParseGapStr(OpenStr, OpenA, LOpenA, ROpenA, OpenB, LOpenB, ROpenB);
	ParseGapStr(ExtStr, ExtA, LExtA, RExtA, ExtB, LExtB, RExtB);
}
/******************************************************************************/
void AlnParams::setNucSubstMx(double Match, double Mismatch) {
	
	unsigned N = unsigned(strlen(Alphabet));

	SubstMx = new MxFloatMatrix(256, 256, 0);
	
	for (unsigned i = 0; i < N; ++i) {
		for (unsigned j = 0; j < N; ++j) {
			float v = float(i == j ? Match : Mismatch);

			Byte ui = (Byte) toupper(Alphabet[i]);
			Byte uj = (Byte) toupper(Alphabet[j]);
			Byte li = (Byte) tolower(ui);
			Byte lj = (Byte) tolower(uj);
			ui = (Byte) toupper(ui);
			uj = (Byte) toupper(uj);

			SubstMx->matrix[ui][uj] = v;
			SubstMx->matrix[uj][ui] = v;

			SubstMx->matrix[ui][lj] = v;
			SubstMx->matrix[uj][li] = v;

			SubstMx->matrix[li][uj] = v;
			SubstMx->matrix[lj][ui] = v;

			SubstMx->matrix[li][lj] = v;
			SubstMx->matrix[lj][li] = v;
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

		SubstMx->matrix[ui][uj] = v;
		SubstMx->matrix[uj][ui] = v;

		SubstMx->matrix[ui][lj] = v;
		SubstMx->matrix[uj][li] = v;

		SubstMx->matrix[li][uj] = v;
		SubstMx->matrix[lj][ui] = v;

		SubstMx->matrix[li][lj] = v;
		SubstMx->matrix[lj][li] = v;
	}
}
/******************************************************************************/
void AlnParams::InitFromCmdLine() {
	Clear();
	
	setNucSubstMx(1.0, -2.0); // opt_match, opt_mismatch

	// Local
	setLocal(-10.0f, -1.0f);
	
	// Global
	OpenA = OpenB = -10.0;
	
	LOpenA = LOpenB = ROpenA = ROpenB = -0.5;
	ExtA = ExtB = -1.0;
	LExtA = LExtB = RExtA = RExtB = -0.5;

	setPenalties("", "");
}
/******************************************************************************/
