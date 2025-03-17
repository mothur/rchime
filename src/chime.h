#ifndef chime_h
#define chime_h

#include "seq.h"
#include "myopts.h"

class ChimeHit2 {

public:

	std::string QLabel;
    std::string ALabel;
    std::string BLabel;
    std::string Q3;
    std::string A3;
    std::string B3;

	//unsigned LY, LN, LA, LD;
	//unsigned RY, RN, RA, RD;
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
		Clear();
	    Options* opt = Options::getInstance();
	    minh = opt->getMinh();
	    mindiv = opt->getMindiv();
	    mindiffs = opt->getMindiffs();
	}

	void Clear() {
		Q3.clear();
		A3.clear();
		B3.clear();
		QLabel.clear();
		ALabel.clear();
		BLabel.clear();

		//LY = LN = LA = LD = UINT_MAX;
		//RY = RN = RA = RD = UINT_MAX;
		ColLo = ColHi = QXLo = QXHi = ColXLo = ColXHi = UINT_MAX;
		CS_LY = CS_LN = CS_LA = CS_RY = CS_RN = CS_RA = UINT_MAX;
		PctIdQT = PctIdQA = PctIdQB = PctIdQM = PctIdAB = -1.0;
		Div = -1.0;
		H = -1.0;
		Score = -1.0;
		AbQ = AbA = AbB = -1.0f;
    }

	bool Accept() const {
		return Score >= minh && Div >= mindiv && CS_LY >= mindiffs && CS_RY >= mindiffs;
	}


private:

    double minh, mindiv, mindiffs;

    bool operator<(const ChimeHit2 &rhs) const {
        if (Score == rhs.Score)
            return Div > rhs.Div;
        return Score > rhs.Score;
    }

};

#endif // chime_h
