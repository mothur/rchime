#ifndef myutils_h
#define myutils_h

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

#include "uchime.h"
#include <random>
#include "mothur.h"

/******************************************************************************/
class Utilities {

public:

    static Utilities* getInstance();

    double Pct(double x, double y);
    unsigned getRandomIndex(unsigned maxIndex);

    void* mymalloc(unsigned n);
    void myfree(void *p, unsigned Bytes);

private:

    static Utilities* _uniqueInstance;
    Utilities( const Utilities& ); // Disable copy constructor
    void operator=( const Utilities& ); // Disable assignment operator

    // set to default options
    Utilities() {}
    ~Utilities() = default;

    double GetMemUseBytes();
    double GetPeakMemUseBytes();

    vector<string> warnings, errors;
    mt19937_64 mersenne_twister_engine;

    double g_PeakMemUseBytes;
    unsigned g_NewCalls;
    unsigned g_FreeCalls;
    double g_InitialMemUseBytes;
    double g_TotalAllocBytes;
    double g_TotalFreeBytes;
    double g_NetBytes;
    double g_MaxNetBytes;
};
/******************************************************************************/




#endif	// myutils_h
