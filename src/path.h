#ifndef path_h
#define path_h

#include "uchime.h"

/******************************************************************************/
struct PathData {

public:
	
	PathData(); 
	~PathData();
	
	void Free();
	void Alloc(unsigned MaxLen);

	char *Start;
	char *Front;
	char *Back;
	unsigned Bytes;

private:

	PathData(PathData &);
	PathData &operator=(PathData &);
};

/******************************************************************************/
const unsigned PathMagic = 0x9A783A16;

struct PathBuffer {
	unsigned Magic;
	char *Buffer;
	unsigned Size;
	bool InUse;
};
/******************************************************************************/
static PathBuffer **g_PathBuffers;
static unsigned g_PathBufferSize;
/******************************************************************************/

#endif // path_h
