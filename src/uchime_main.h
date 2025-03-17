#ifndef UCHIME_MAIN_H
#define UCHIME_MAIN_H

/*
 *  uchime_main.h
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "myutils.h"
#include "chime.h"
#include "seqdb.h"
#include "searchchime.h"
#include "uchime_output_converter.h"
#include "mothur-r.h"

//const vector<float> *g_SortVecFloat;
//bool g_UchimeDeNovo = false;

/* This class is the access point to the uchime source */
/******************************************************************************/
class UchimeMain {

public:

    UchimeMain();   
    ~UchimeMain();  

    // returns list of dataframes with results
    Rcpp::List runUchime(vector<string>& names, vector<string>& seqs,
                         vector<string>& refNames, vector<string>& refSeqs,
                         vector<int>& abunds, bool chimeAlns);

private:

   // create a dataframe to pass back to R
   UchimeOutputConverter* dataUchime;
   SeqDB* reference;
   SeqDB* data;
   Utils* util;

   bool uchimeDeNovo;
};
/******************************************************************************/
#endif
