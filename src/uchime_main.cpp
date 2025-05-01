
/*
 *  uchime_main.cpp
 *
 *  Modifications by Sarah Westcott on 2/20/25.
 *
 */

#include "uchime_main.h"

/******************************************************************************/
// this is the entry point to uchime source code from chimera_findur
vector<ChimeHit2> UchimeMain::runUchime(vector<string> names,
                                 vector<string> seqs,
                                 vector<string> refNames,
                                 vector<string> refSeqs,
                                 vector<int> abunds,
								 set<string>& namesOfChimeras,
								 Options opts, bool silent, SEXP* bar) {

	// loading SeqDB with data from R, sort descending order by abundace
	SeqDB data(names, seqs, abunds, true, false);

	// are we running with a reference, or denovo
	uchimeDeNovo = (refNames.size() == 0);

	SeqDB* reference = nullptr;
	if (uchimeDeNovo) {
		reference = new SeqDB(true, true);
	}else{
		vector<int> refAbunds(refNames.size(), 1);
		reference = new SeqDB(refNames, refSeqs, refAbunds, false, true);
	}

	vector<ChimeHit2> Hits;
	unsigned numQuerySeqs = data.getSeqCount();
	SearchChime search(opts);
	for (unsigned i = 0; i < numQuerySeqs; ++i) {

		// main thread updates progress
		if (!silent && (i % 500 == 0)) {
			cli_progress_set(*bar, i+1);
		}

		SeqData queryData = data.getSeqData(i);

		ChimeHit2 Hit;
		bool chimeric = search.searchChime(reference, queryData, Hit);
		if (!chimeric) {
			// add more abundant query sequences as potential parents
			if (uchimeDeNovo) {
				reference->addSeq(queryData.getName(), queryData.getSeq(), queryData.getAbund());
			}
		}

		Hits.push_back(Hit);
		if (chimeric) {  namesOfChimeras.insert(queryData.getName()); }
	}

	delete reference;
	return Hits;
}
/******************************************************************************/
