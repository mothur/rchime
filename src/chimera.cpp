
/*
 *  chimera.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "chimera.h"

Chimera::Chimera(bool derep, int proc, bool si) {
    dereplicate = derep;
    processors = proc;
    silent = si;
}
