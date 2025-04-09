#ifndef chime_h
#define chime_h

#include "seq.h"
#include "options.h"

struct ChimeHit2 {

	string QLabel;
    string ALabel;
    string BLabel;
    string Q3;
    string A3;
    string B3;

	double PctIdQT, PctIdQA, PctIdQB, PctIdQM, PctIdAB;

	unsigned ColLo;
	unsigned ColXLo;
	unsigned ColXHi;
	unsigned ColHi;
	unsigned QXLo;
	unsigned QXHi;

	double Div;
	double Score;
	double H;

	unsigned CS_LY, CS_LN, CS_LA, CS_RY, CS_RN, CS_RA;

	float AbQ;
	float AbA;
	float AbB;

	ChimeHit2() {
		Q3.clear();
		A3.clear();
		B3.clear();
		QLabel.clear();
		ALabel.clear();
		BLabel.clear();

		
		ColLo = ColHi = QXLo = QXHi = ColXLo = ColXHi = UINT_MAX;
		CS_LY = CS_LN = CS_LA = CS_RY = CS_RN = CS_RA = UINT_MAX;
		PctIdQT = PctIdQA = PctIdQB = PctIdQM = PctIdAB = -1.0;
		Div = -1.0;
		H = -1.0;
		Score = -1.0;
		AbQ = AbA = AbB = -1.0f;

	    Options* opt = Options::getInstance();
	    minh = opt->getMinh();
	    mindiv = opt->getMindiv();
	    mindiffs = opt->getMindiffs();
	}

	bool Accept() const {
		return Score >= minh && Div >= mindiv && CS_LY >= mindiffs && CS_RY >= mindiffs;
	}

private:

    double minh, mindiv, mindiffs;

    bool operator<(const ChimeHit2 &rhs) const {
        if (Score == rhs.Score) {
            return Div > rhs.Div;
		}
        return Score > rhs.Score;
    }

};

#endif // chime_h
