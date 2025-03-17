#ifndef CHIMERAUCHIME_H
#define CHIMERAUCHIME_H

/*
 *  chimera_uchime.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

/* This class is inherits from Chimera */

#include "mothur.h"
#include "mothur-r.h"
#include "chimera.h"


/******************************************************************************/
class ChimeraUchime : public Chimera {

public:

    // dereplicate, processors, silent
    ChimeraUchime(bool derep, int proc, bool si);

    Rcpp::List removeChimeras(Rcpp::Environment& dataset);
    Rcpp::List removeChimeras(Rcpp::Environment& dataset,
                                      Rcpp::Environment& reference);

private:


};
/******************************************************************************/
#endif
