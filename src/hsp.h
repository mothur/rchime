#ifndef hsp_h
#define hsp_h	1

#include "uchime.h"

struct HSPData
	{
	unsigned Loi;
	unsigned Loj;
	unsigned Leni;
	unsigned Lenj;
	float Score;
	unsigned User;

	unsigned GetLength() const  { return Leni;           }
	unsigned GetHii() const     { return Loi + Leni - 1; }
	unsigned GetHij() const     { return Loj + Lenj - 1; }

	bool LeftA() const          { return Loi == 0;       }
	bool LeftB() const          { return Loj == 0;       }

	bool RightA(unsigned LA) const { return Loi + Leni == LA; }
	bool RightB(unsigned LB) const { return Loj + Lenj == LB; }

	unsigned GetIdCount(const Byte *A, const Byte *B) const {
		unsigned Count = 0;
		unsigned K = GetLength();

		for (unsigned k = 0; k < K; ++k) {
			Byte a = A[Loi+k];
			Byte b = B[Loj+k];
			if (toupper(a) == toupper(b)) {
				Count++;
			}
		}
		return Count;
	}

	double OverlapFract(const HSPData &HSP) const {
		if (Leni == 0 || Lenj == 0) {
			return 0.0;
		}
		unsigned MaxLoi = max(Loi, HSP.Loi);
		unsigned MaxLoj = max(Loj, HSP.Loj);
		unsigned MinHii = min(GetHii(), HSP.GetHii());
		unsigned MinHij = min(GetHij(), HSP.GetHij());

		unsigned Ovi = (MinHii < MaxLoi) ? 0 : MinHii - MaxLoi;
		unsigned Ovj = (MinHij < MaxLoj) ? 0 : MinHij - MaxLoj;

		return double(Ovi*Ovj)/double(Leni*Lenj);
	}

	bool operator<(const HSPData &rhs) const {
		return Loi < rhs.Loi;
	}
};

// Bendpoint
struct BPData
	{
	unsigned Pos;
	bool IsLo;
	unsigned Index;
};

#endif // hsp_h
