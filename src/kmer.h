#ifndef KMER_HPP
#define KMER_HPP

/*
 *  kmer.hpp
 *
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"

/**************************************************************************************************/

class Kmer {

public:
	Kmer(int);
    ~Kmer() = default;

	int getKmerNumber(string, int);

private:

	int kmerSize;
	int maxKmer;
};

/**************************************************************************************************/


#endif
