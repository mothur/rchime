
/*
 *  The Utils class is based on myutils.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 2/24/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will handle warnings, errors, memory management, and minor utility
 * functions.
 *
 */

#include "myutils.h"

/******************************************************/
Utils* Utils::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new Utils();
	}
	return _uniqueInstance;
}
/******************************************************************************/
unsigned Utils::getRandomIndex(unsigned maxIndex){
    if (maxIndex == 0) {  return 0; }

    uniform_int_distribution<unsigned> dis(0, maxIndex);

    unsigned random = dis(mersenne_twister_engine);

    return (random);
}
/******************************************************************************/
void* Utils::mymalloc(unsigned Bytes) {
    ++g_NewCalls;
    if (g_InitialMemUseBytes == 0)
        g_InitialMemUseBytes = GetMemUseBytes();

    g_TotalAllocBytes += Bytes;
    g_NetBytes += Bytes;
    if (g_NetBytes > g_MaxNetBytes)
    {
        if (g_NetBytes > g_MaxNetBytes + 10000000)
            GetMemUseBytes();//to force update of peak
        g_MaxNetBytes = g_NetBytes;
    }
    void *p = malloc(Bytes);

    // // if we were unable to allocate the memory, send error message
    // if (p == nullptr) {
    //     double b = GetMemUseBytes();

    //     string errorMessage = "Out of memory, attempted allocating " +
    //         toString(Bytes) + ". Current usage " + toString(b) + ".";

    //     Die(errorMessage);
    // }
    return p;
}
/******************************************************************************/
void Utils::myfree(void *p, unsigned Bytes) {
    ++g_FreeCalls;
    g_TotalFreeBytes += Bytes;
    g_NetBytes -= Bytes;

    if (p != nullptr) {
        free(p);
    }
}
/******************************************************************************/
double Utils::Pct(double x, double y) {
    if (y == 0.0f) {
        return 0.0f;
    }
    return (x*100.0f)/y;
}
/******************************************************************************/
double Utils::GetPeakMemUseBytes() {
    return g_PeakMemUseBytes;
}
/******************************************************************************/
#ifdef _MSC_VER

double Utils::GetMemUseBytes() {
	HANDLE hProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS PMC;
	BOOL bOk = GetProcessMemoryInfo(hProc, &PMC, sizeof(PMC));
	if (!bOk)
		return 1000000;
	double Bytes = (double) PMC.WorkingSetSize;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
}

#elif	linux || __linux__

double Utils::GetMemUseBytes() {
	static char statm[64];
	static int PageSize = 1;
	if (0 == statm[0])
		{
		PageSize = sysconf(_SC_PAGESIZE);
		pid_t pid = getpid();
		sprintf(statm, "/proc/%d/statm", (int) pid);
		}

	int fd = open(statm, O_RDONLY);
	if (-1 == fd)
		return 1000000;
	char Buffer[64];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		return 1000000;

	Buffer[n] = 0;
	double Pages = atof(Buffer);

	double Bytes = Pages*PageSize;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
}

#elif defined(__MACH__)

#define DEFAULT_MEM_USE	100000000.0

double Utils::GetMemUseBytes() {
	task_t mytask = mach_task_self();
	struct task_basic_info ti;
	memset((void *) &ti, 0, sizeof(ti));
	mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;
	kern_return_t ok = task_info(mytask, TASK_BASIC_INFO, (task_info_t) &ti, &count);
	if (ok == KERN_INVALID_ARGUMENT)
		return DEFAULT_MEM_USE;

	if (ok != KERN_SUCCESS)
		return DEFAULT_MEM_USE;

	double Bytes = (double ) ti.resident_size;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
}

#else

double Utils::GetMemUseBytes() {
	return 0;
}

#endif

/******************************************************************************/
