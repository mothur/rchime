#ifndef UCHIME_h
#define UCHIME_h

#include <stdio.h>
#include <sys/types.h>
#include <string>
#include <string.h>
#include <memory.h>
#include <vector>
#include <math.h>
#include <stdarg.h>
#include <cstdlib>
#include <climits>
#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>
#include <set>
#include <map>
#include <signal.h>
#include <float.h>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <float.h>
#include <list>
#include <limits.h>


#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#else
#include <crtdbg.h>
#include <process.h>
#include <windows.h>
#include <psapi.h>
#include <io.h>
#include <inttypes.h>
#endif

#ifdef __MACH__

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/socket.h>
#include <sys/gmon.h>
#include <mach/vm_param.h>
#include <netinet/in.h>
#include <netinet/icmp6.h>
#include <sys/vmmeter.h>
#include <sys/proc.h>
#include <mach/vm_statistics.h>
#include <mach/task_info.h>
#include <mach/task.h>
#include <mach/mach_init.h>

#endif

using namespace std;

/******************************************************************************/
// inline functions
// pom=plus or minus, tof=true or false
static inline char pom(bool Plus)	{ return Plus ? '+' : '-'; }
static inline char tof(bool x)		{ return x ? 'T' : 'F';	}
static inline char yon(bool x)		{ return x ? 'Y' : 'N';	}

// Are two floats equal to within epsilon?
const double epsilon = 0.01;
const unsigned wordCount = 65536;

struct orderAbundance {
    int abund;
    unsigned index;

    orderAbundance() { abund = 0; index = 0; }
    orderAbundance(int a, unsigned i) : abund(a), index(i) {}
    ~orderAbundance() {}
};

struct orderFloatAbundance {
    float abund;
    unsigned index;

    orderFloatAbundance() { abund = 0; index = 0; }
    orderFloatAbundance(float a, unsigned i) : abund(a), index(i) {}
    ~orderFloatAbundance() {}
};

//sorts highest to lowest
static inline bool compareAbundance(orderAbundance left, orderAbundance right){
	return (left.abund > right.abund);	
} 

static inline bool compareFloatAbundance(orderFloatAbundance left, orderFloatAbundance right){
	return (left.abund > right.abund);	
} 

template<typename T>
void applyOrder(T& x, const std::vector<unsigned>& order) {
    T copy = x;

    for (int i = 0; i < order.size(); i++) {
        x[i] = copy[order[i]];
    }
}

static inline bool isacgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

static bool inline isgap(char c) {
    return c == '-' || c == '.';
}
/******************************************************************************/
// definitions

#ifdef _MSC_VER
#include <crtdbg.h>
#pragma warning(disable: 4996)	// deprecated functions
#define _CRT_SECURE_NO_DEPRECATE	1
#endif

typedef unsigned char Byte;
typedef unsigned short uint16;
typedef unsigned uint32;
typedef int int32;
typedef double float32;
typedef signed char int8;
typedef unsigned char uint8;

#ifdef _MSC_VER

typedef __int64 int64;
typedef unsigned __int64 uint64;

#elif defined(__x86_64__)

typedef long int64;
typedef unsigned long uint64;

#else

typedef long long int64;
typedef unsigned long long uint64;

#endif

#ifdef _MSC_VER
#define off_t	__int64
#endif

#define SIZE(c)	unsigned((c).size())

/******************************************************************************/

#endif
