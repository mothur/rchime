#ifndef path_h
#define path_h

#include "uchime.h"
#include "myutils.h"

/******************************************************************************/
struct PathData {

private:
	PathData(PathData &);
	PathData &operator=(PathData &);

public:
	char *Start;
	char *Front;
	char *Back;
	unsigned Bytes;

public:
	PathData() {
		Clear(true);
	}
	~PathData() {
		Free();
	 }
	void Free();
	void Alloc(unsigned MaxLen);
	void Clear(bool ctor = false) {
		Start = 0;
		if (ctor) {
			Front = 0;
			Back = 0;
			Bytes = 0;
		}else {
			Free();
		}
	}
};

/******************************************************************************/
const unsigned PathMagic = 0x9A783A16;

struct PathBuffer {
	unsigned Magic;
	char *Buffer;
	unsigned Size;
	bool InUse;
};

static PathBuffer **g_PathBuffers;
static unsigned g_PathBufferSize;
/******************************************************************************/

#endif // path_h
