#ifndef alpha_h
#define alpha_h

#include <limits.h>
#include <string>

using namespace std;

const unsigned INVALID_LETTER = 0;
const unsigned char INVALID_CHAR = '?';

extern unsigned g_CharToLetterNucleo[];
extern bool g_IsNucleoChar[];
extern bool g_IsACGTU[];

#endif // alpha_h
