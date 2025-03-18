
#include "path.h"

/******************************************************************************/
static char *AllocBuffer(unsigned Size){
	if (Size == 0) {
		return 0;
	}
	Utilities* util = Utilities::getInstance();
	
	// Is a free buffer that is big enough?
	for (unsigned i = 0; i < g_PathBufferSize; ++i) {
		PathBuffer *PB = g_PathBuffers[i];
		
		if (!PB->InUse) {
			if (PB->Size >= Size) {
				PB->InUse = true;
				return PB->Buffer;
			}
			if (PB->Buffer == 0) {
				unsigned Size2 = Size + 1024;
				PB->Buffer = static_cast<char*> (util->mymalloc(sizeof(char) * Size2));
				PB->Size = Size2;
				PB->InUse = true;
				return PB->Buffer;
			}
		}
	}

	// No available buffer, must expand g_PathBuffers[]
	unsigned NewPathBufferSize = g_PathBufferSize + 1024;
	PathBuffer **NewPathBuffers = static_cast<PathBuffer**> (util->mymalloc(sizeof(PathBuffer *) * NewPathBufferSize));
	
	for (unsigned i = 0; i < g_PathBufferSize; ++i) {
		NewPathBuffers[i] = g_PathBuffers[i];
	}

	for (unsigned i = g_PathBufferSize; i < NewPathBufferSize; ++i) {
		PathBuffer *PB = static_cast<PathBuffer*> (util->mymalloc(sizeof(PathBuffer) * 1)); 
		PB->Magic = PathMagic;
		PB->Buffer = 0;
		PB->Size = 0;
		PB->InUse = false;
		NewPathBuffers[i] = PB;
	}

	PathBuffer *PB = NewPathBuffers[g_PathBufferSize];

	util->myfree(g_PathBuffers, g_PathBufferSize);
	g_PathBuffers = NewPathBuffers;
	g_PathBufferSize = NewPathBufferSize;

	unsigned Size2 = Size + 1024;
	PB->Buffer = static_cast<char*> (util->mymalloc(sizeof(char) * Size2)); 
	PB->Size = Size2;
	PB->InUse = true;
	return PB->Buffer;
}
/******************************************************************************/
static void FreeBuffer(char *Buffer) {
	if (Buffer == 0) {
		return;
	}

	for (unsigned i = 0; i < g_PathBufferSize; ++i){
		PathBuffer *PB = g_PathBuffers[i];
		if (PB->Buffer == Buffer){
			PB->InUse = false;
			return;
		}
	}
}
/******************************************************************************/
void PathData::Alloc(unsigned MaxLen) {
	if (MaxLen < Bytes) {
		return;
	}
	if (Bytes > 0) {
		FreeBuffer(Front);
	}

	Bytes = MaxLen + 1;
	Front = AllocBuffer(Bytes);
	Back = Front + Bytes - 1;
	Start = 0;
}
/******************************************************************************/
void PathData::Free() {
	FreeBuffer(Front);
	Front = 0;
	Start = 0;
	Back = 0;
}
/******************************************************************************/
