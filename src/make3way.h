#ifndef MAKE3WAY_H
#define MAKE3WAY_H

/*
 *  The Make3way class is based on make3way.cpp from Robert Edgar.
 *
 *  Modified by Sarah Westcott on 2/26/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime.h"
#include "myutils.h"
#include "chime.h"

class Make3way {
    public:

        Make3way() = default;
        ~Make3way() = default;


        void Make3Way(const SeqData &QSD, const SeqData &ASD, const SeqData &B_SD,
  						const string &PathQA, const string &PathQB, string &Q3, 
						string &A3, string &B3);

};

#endif
