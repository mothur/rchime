#include "myutils.h"
#include "mx.h"
#include "seqdb.h"
#include "seq.h"

char ProbToChar(float p);

/******************************************************************************/
list<MxBase *> *MxBase::m_Matrices = 0;
unsigned MxBase::m_AllocCount;
unsigned MxBase::m_ZeroAllocCount;
unsigned MxBase::m_GrowAllocCount;
double MxBase::m_TotalBytes;
double MxBase::m_MaxBytes;
/******************************************************************************/
void MxBase::OnCtor(MxBase *Mx) {
	if (m_Matrices == 0) {
		m_Matrices = new list<MxBase *>;
	}
	
	m_Matrices->push_front(Mx);
}
/******************************************************************************/
void MxBase::OnDtor(MxBase *Mx) {
	if (m_Matrices == 0) {
		return;
	}

	for (list<MxBase*>::iterator p = m_Matrices->begin();
	  p != m_Matrices->end(); ++p) {
		if (*p == Mx) {
			m_Matrices->erase(p);
			if (m_Matrices->empty()) {
				delete m_Matrices;
			}
			return;
		}
	}
}
/******************************************************************************/
void MxBase::Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
  const SeqDB *DB, unsigned IdA, unsigned IdB) {
	Alloc(Name, RowCount, ColCount, DB, IdA, IdB, 0, 0);
}
/******************************************************************************/
void MxBase::Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
  const SeqData *SA, const SeqData *SB) {
	Alloc(Name, RowCount, ColCount, 0, UINT_MAX, UINT_MAX, SA, SB);
}
/******************************************************************************/
void MxBase::Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
  const SeqDB *DB, unsigned IdA, unsigned IdB, const SeqData *SA, const SeqData *SB) {

	++m_AllocCount;
	if (m_AllocatedRowCount == 0) {
		++m_ZeroAllocCount;
	}

	if (RowCount > m_AllocatedRowCount || ColCount > m_AllocatedColCount) {
		if (m_AllocatedRowCount > 0) {
			++m_GrowAllocCount;
		}

		m_TotalBytes -= GetBytes();
		FreeData();

		unsigned N = max(RowCount + 16, m_AllocatedRowCount);
		unsigned M = max(ColCount + 16, m_AllocatedColCount);
		N = max(N, M);

		AllocData(N, N);
		
		m_TotalBytes += GetBytes();
		if (m_TotalBytes > m_MaxBytes) {
			m_MaxBytes = m_TotalBytes;
		}
	}

	unsigned n = sizeof(m_Name)-1;
	strncpy(m_Name, Name, n);
	m_Name[n] = 0;
	m_RowCount = RowCount;
	m_ColCount = ColCount;
	m_SeqDB = DB;
	m_IdA = IdA;
	m_IdB = IdB;
	m_SA = SA;
	m_SB = SB;
}
/******************************************************************************/
//static unsigned g_MatrixFileCount;
