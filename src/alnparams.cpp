
#include "alnparams.h"
#include "hsp.h"

/******************************************************************************/
AlnParams::AlnParams() {   
	opt = Options::getInstance();
	util = Utilities::getInstance();
	InitFromCmdLine(true);
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
	Nucleo = false;
	NucleoSet = false;
}
/******************************************************************************/
void AlnParams::setLocal(float Open, float Ext) {
	LocalOpen = Open;
	LocalExt = Ext;
}
/******************************************************************************/
void AlnParams::Init4(float Open, float Ext, float TermOpen, float TermExt) {
	SubstMx = g_SubstMx;
	OpenA = OpenB = Open;
	LOpenA = LOpenB = ROpenA = ROpenB = TermOpen;
	ExtA = ExtB = Ext;
	LExtA = LExtB = RExtA = RExtB = TermExt;
}
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
static void ParseGapStr(const string &s, float &QI, float &QL, float &QR,
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
void AlnParams::SetMxFromCmdLine(bool IsNucleo) {
	if (IsNucleo) {
		setNucSubstMx(1.0, -2.0); // opt_match, opt_mismatch
	}
	SubstMx = g_SubstMx;
}
/******************************************************************************/
void AlnParams::setNucSubstMx(double Match, double Mismatch) {
	
	unsigned N = unsigned(strlen(Alphabet));

	g_SubstMxf.Alloc("NUCMX", 256, 256);
	strcpy(g_SubstMxf.m_Alpha, "ACGT");
	g_SubstMxf.Init(0);
	g_SubstMx = g_SubstMxf.GetData();

	for (unsigned i = 0; i < N; ++i) {
		for (unsigned j = 0; j < N; ++j) {
			float v = float(i == j ? Match : Mismatch);

			Byte ui = (Byte) toupper(Alphabet[i]);
			Byte uj = (Byte) toupper(Alphabet[j]);
			Byte li = (Byte) tolower(ui);
			Byte lj = (Byte) tolower(uj);
			ui = (Byte) toupper(ui);
			uj = (Byte) toupper(uj);

			g_SubstMx[ui][uj] = v;
			g_SubstMx[uj][ui] = v;

			g_SubstMx[ui][lj] = v;
			g_SubstMx[uj][li] = v;

			g_SubstMx[li][uj] = v;
			g_SubstMx[lj][ui] = v;

			g_SubstMx[li][lj] = v;
			g_SubstMx[lj][li] = v;
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

		g_SubstMx[ui][uj] = v;
		g_SubstMx[uj][ui] = v;

		g_SubstMx[ui][lj] = v;
		g_SubstMx[uj][li] = v;

		g_SubstMx[li][uj] = v;
		g_SubstMx[lj][ui] = v;

		g_SubstMx[li][lj] = v;
		g_SubstMx[lj][li] = v;
	}
}
/******************************************************************************/
void AlnParams::InitFromCmdLine(bool IsNucleo) {
	Clear();
	Nucleo = IsNucleo;
	NucleoSet = true;

	SetMxFromCmdLine(IsNucleo);

	// Local
	setLocal(-10.0f, -1.0f);
	
	// Global
	if (IsNucleo) {
		Init4(-10.0, -1.0, -0.5, -0.5);
	}//else {
		//Init4(-17.0, -1.0, -0.5, -0.5);
	//}
	setPenalties(0, 0);
}
/******************************************************************************/
