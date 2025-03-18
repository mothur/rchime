#ifndef mx_h
#define mx_h

#include "uchime.h"
#include "myutils.h"
#include "seqdb.h"

const int OPT_LOG = 0x01;
const int OPT_EXP = 0x02;
const int OPT_ZERO_BASED = 0x04;
const float MINUS_INFINITY = -9e9f;
const float UNINIT = -8e8f;

/******************************************************************************/
struct MxBase {
private:
	MxBase(const MxBase &rhs);
	MxBase &operator=(const MxBase &rhs);

public:
	char m_Name[32];
	char m_Alpha[32];
	char m_Alpha2[32];
	unsigned m_RowCount;
	unsigned m_ColCount;
	unsigned m_AllocatedRowCount;
	unsigned m_AllocatedColCount;
	const SeqDB *m_SeqDB;
	unsigned m_IdA;
	unsigned m_IdB;
	const SeqData *m_SA;
	const SeqData *m_SB;

	static list<MxBase *> *m_Matrices;

	static unsigned m_AllocCount;
	static unsigned m_ZeroAllocCount;
	static unsigned m_GrowAllocCount;
	static double m_TotalBytes;
	static double m_MaxBytes;

	static void OnCtor(MxBase *Mx);
	static void OnDtor(MxBase *Mx);

	MxBase() {
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_IdA = UINT_MAX;
		m_IdB = UINT_MAX;
		m_SeqDB = 0;
		OnCtor(this);
	}
	virtual ~MxBase() {
		OnDtor(this);
	}

	virtual unsigned GetTypeSize() const = 0;
	virtual unsigned GetBytes() const = 0;

	void Clear() {
		FreeData();
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_IdA = UINT_MAX;
		m_IdB = UINT_MAX;
		m_SA = 0;
		m_SB = 0;
	}

	bool Empty() const {
		return m_RowCount == 0;
	}

	virtual void AllocData(unsigned RowCount, unsigned ColCount) = 0;
	virtual void FreeData() = 0;

	void SetAlpha(const char *Alpha) {
		unsigned n = sizeof(m_Alpha);
		strncpy(m_Alpha, Alpha, n);
		m_Alpha[n] = 0;
	}

	void Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
	  const SeqDB *DB, unsigned IdA, unsigned IdB,
	  const SeqData *SA, const SeqData *SB);

	void Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
	  const SeqDB *DB = 0, unsigned IdA = UINT_MAX, unsigned IdB = UINT_MAX);

	void Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
	  const SeqData *SA, const SeqData *SB);
};
/******************************************************************************/
template<class T> struct Mx : public MxBase {
	
// Disable unimplemented stuff
private:
	Mx(Mx &rhs);
	Mx &operator=(Mx &rhs);
	
    Utilities* util;

public:
	T **m_Data;

	Mx()  { m_Data = 0; util = Utilities::getInstance(); }
	~Mx() { FreeData(); }

	virtual void AllocData(unsigned RowCount, unsigned ColCount) {
		m_Data = static_cast<T**> (util->mymalloc(sizeof(T*) * RowCount));
		for (unsigned i = 0; i < RowCount; ++i)
			// m_Data[i] = myalloc<T>(ColCount);
			m_Data[i] = static_cast<T*> (util->mymalloc(sizeof(T) * ColCount));

		m_AllocatedRowCount = RowCount;
		m_AllocatedColCount = ColCount;
	}

	virtual void FreeData() {
		for (unsigned i = 0; i < m_AllocatedRowCount; ++i) {
			util->myfree(m_Data[i], m_AllocatedColCount);
        }
		util->myfree(m_Data, m_AllocatedRowCount);

		m_Data = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
	}

	T **GetData() {
		return (T **) m_Data;
	}

	unsigned GetTypeSize() const {
		return sizeof(T);
	}

	virtual unsigned GetBytes() const {
		return m_AllocatedRowCount*m_AllocatedColCount*GetTypeSize() +
		  m_AllocatedRowCount*sizeof(T *);
	}

	const T *const *const GetData() const {
		return (const T *const *) m_Data;
	}

	void Init(T v) {
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = v;
	}
};
/******************************************************************************/

#endif // mx_h
