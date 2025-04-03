
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
#include <__config>
#include <cli/progress.h>

/******************************************************************************/
Options* Options::_uniqueInstance = 0;
Utilities* Utilities::_uniqueInstance = 0;
/******************************************************************************/
UchimeMain::UchimeMain()    { 
	util = Utilities::getInstance();
}
/******************************************************************************/
UchimeMain::~UchimeMain()   {}

/******************************************************************************/
// this is the entry point to uchime source code from chimera_findur
vector<ChimeHit2> UchimeMain::runUchime(vector<string> names,
                                 vector<string> seqs,
                                 vector<string> refNames,
                                 vector<string> refSeqs,
                                 vector<int> abunds,
								 set<string>& namesOfChimeras)
	{

    vector<ChimeHit2> Hits;

	// loading SeqDB with data from R, sort descending order by abundace
	SeqDB data(names, seqs, abunds, true);

	// are we running with a reference, or denovo
	uchimeDeNovo = (refNames.size() == 0);

	SeqDB* reference = nullptr;
	if (uchimeDeNovo) {
		reference = new SeqDB();
	}else{
		vector<int> refAbunds(refNames.size(), 1);
		reference = new SeqDB(refNames, refSeqs, refAbunds, false);
	}

	unsigned numQuerySeqs = data.getSeqCount();
	SearchChime search;
	for (unsigned i = 0; i < numQuerySeqs; ++i) {

		SeqData queryData = data.getSeqData(i);
if (queryData.getName() == "M00967_43_000000000-A3JHG_1_1114_16114_2126" ) {
    cout << "here;";
}
		ChimeHit2 Hit;
		bool chimeric = search.searchChime(reference, queryData, Hit);
		if (!chimeric) {
			// add more abundant query sequences as potential parents
			if (uchimeDeNovo) {
				reference->addSeq(queryData.getName(), queryData.getSeq(), queryData.getAbund(), false);
			}
		}

		Hits.push_back(Hit);
		if (chimeric) {  namesOfChimeras.insert(queryData.getName()); }
	}

	delete reference;
	return Hits;
}
/******************************************************************************/
