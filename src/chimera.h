#ifndef CHIMERA_H
#define CHIMERA_H

/*
 *  chimera.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

/* This class is a parent to ChimeraUchime, ChimeraVsearch */

#include "mothur.h"
#include "mothur-r.h"
#include "rcpp_utils.h"
#include "myopts.h"

/******************************************************************************/
class Chimera {

public:

    // dereplicate, processors, silent
    Chimera(bool derep, int proc, bool si, bool hg, Rcpp::List options);

    virtual ~Chimera(){};

    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset) = 0;
    virtual Rcpp::List removeChimeras(Rcpp::Environment& dataset,
                                      Rcpp::Environment& reference) = 0;

protected:

    int processors;
    bool silent, dereplicate, hasGroupData;
    Options* opts;

private:


};
/******************************************************************************/
#endif
