
/*
 *  uchime_main.cpp
 *
 *  Modifications by Sarah Westcott on 2/20/25.
 *
 */

#include "uchime_main.h"
#include "myopts.h"
#include "myutils.h"
#include "seqdb.h"
#include <Rcpp.h>
#include <RcppThread.h>
#include <cli/progress.h>


/******************************************************************************/
Options* Options::_uniqueInstance = 0;
Utilities* Utilities::_uniqueInstance = 0;
/******************************************************************************/
UchimeMain::UchimeMain()    { 
	util = Utilities::getInstance();
	reference = nullptr;
	data = nullptr;
}
/******************************************************************************/
UchimeMain::~UchimeMain()   {
	 if (reference != nullptr) { delete reference; }
	 if (data != nullptr) { delete data; }
}
/******************************************************************************/
// this is the entry point to uchime source code from chimera_findur
vector<ChimeHit2> UchimeMain::runUchime(vector<string>& names,
                                 vector<string>& seqs,
                                 vector<string>& refNames,
                                 vector<string>& refSeqs,
                                 vector<int>& abunds,
								 set<string>& namesOfChimeras)
	{

    vector<ChimeHit2> Hits;

	// loading SeqDB with data from R, sort descending order by abundace
	data = new SeqDB(names, seqs, abunds, true);

	if (!data->isNucleo()) {
	    return Hits;
	}

	// are we running with a reference, or denovo
	uchimeDeNovo = (refNames.size() == 0);

	if (uchimeDeNovo) {
		reference = new SeqDB();
	}else{
		vector<int> refAbunds(refNames.size(), 1);
		reference = new SeqDB(refNames, refSeqs, refAbunds, false);
		if (!reference->isNucleo()) {
		    return Hits;
		}
	}

	unsigned numQuerySeqs = data->getSeqCount();
	SearchChime search;
	for (unsigned i = 0; i < numQuerySeqs; ++i) {

		SeqData queryData = data->getSeqData(i);

		ChimeHit2 Hit;
		bool Found = search.searchChime(reference, queryData, Hit);
		if (!Found) {
			// add more abundant query sequences as potential parents
			if (uchimeDeNovo) {
				reference->addSeq(queryData.getName(), queryData.getSeq(), queryData.getAbund());
			}
		}

		Hits.push_back(Hit);
		if (Hit.Accept()) {  
			namesOfChimeras.insert(queryData.getName());
		}
	}

	return Hits;
}
/******************************************************************************/
