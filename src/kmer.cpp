/*
 *  kmer.cpp
 *
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "kmer.h"

/**************************************************************************************************/

Kmer::Kmer(int size) : kmerSize(size) {	//	The constructor sets the size of the kmer

	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
																		//	No reason to waste the time of recalculating
	maxKmer = power4s[kmerSize]+1;// (int)pow(4.,k)+1;					//	powers of 4 everytime through.  We need an
																		//	extra kmer if we get a non-ATGCU base
}
/**************************************************************************************************/

int Kmer::getKmerNumber(string sequence, int index){

//	Here we convert a kmer to a number between 0 and maxKmer.  For example, AAAA would equal 0 and TTTT would equal 255.
//	If there's an N in the kmer, it is set to 256 (if we are looking at 4mers).  The largest we can look at are 8mers,
//	this could be easily increased.
//
//	Example:   ATGCAT (kSize = 6)
//		  i:   012345
//
//		Score	= 0*4^(6-0-1) + 3*4^(6-1-1) + 2*4^(6-2-1) + 1*4^(6-3-1) + 0*4^(6-4-1) + 3*4^(6-5-1)
//				= 0*4^5	+	3*4^4	+	2*4^3	+	1*4^2	+	0*4^1	+	3*4^0
//				= 0 + 3*256 + 2*64 + 1*16 + 0*4 + 3*1
//				= 0 + 768 + 128 + 16 + 0 + 3
//				= 915

	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };

	int kmer = 0;
	for(int i=0;i<kmerSize;i++){
		if(toupper(sequence[i+index]) == 'A')		{	kmer += (0 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'C')	{	kmer += (1 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'G')	{	kmer += (2 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'U')	{	kmer += (3 * power4s[kmerSize-i-1]);	}
		else if(toupper(sequence[i+index]) == 'T')	{	kmer += (3 * power4s[kmerSize-i-1]);	}
		// invalid letter to match original uchime word calc
		else 	                                    {	kmer += (0 * power4s[kmerSize-i-1]);	}
	}
	return kmer;
}
/**************************************************************************************************/

