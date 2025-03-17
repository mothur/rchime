
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
Utils* Utils::_uniqueInstance = 0;
AlnParams* AlnParams::_uniqueInstance = 0;
/******************************************************************************/
UchimeMain::UchimeMain()    { 
	dataUchime = nullptr; 
	util = Utils::getInstance();
	reference = nullptr;
	data = nullptr;
}
/******************************************************************************/
UchimeMain::~UchimeMain()   {
	 if (dataUchime != nullptr) { delete dataUchime; }
	 if (reference != nullptr) { delete reference; }
	 if (data != nullptr) { delete data; }
}
/******************************************************************************/
// this is the entry point to uchime source code from chimera_findur
Rcpp::List UchimeMain::runUchime(vector<string>& names,
                                 vector<string>& seqs,
                                 vector<string>& refNames,
                                 vector<string>& refSeqs,
                                 vector<int>& abunds, bool chimeAlns)
	{

    Rcpp::List results = Rcpp::List::create();
    vector<string> resultsNames;

    Options* options = Options::getInstance();

	// loading SeqDB with data from R, sort descending order by abundace
	data = new SeqDB(names, seqs, abunds, true);

	if (!data->isNucleo()) {
	    return results;
	}

	// replace file output with object - uchimeout file
	dataUchime = new UchimeOutputConverter();

	if (options->getChimealns()) {
	    // 	g_fUChimeAlns = CreateStdioFile(opt_uchimealns);
	    // WriteChimeHitX creates the file from the Hits
	    // resultsNames.push_back("chimealns");
	}

	// are we running with a reference, or denovo
	uchimeDeNovo = (refNames.size() == 0);

	if (uchimeDeNovo) {
		reference = new SeqDB();
	}else{
		vector<int> refAbunds(refNames.size(), 1);
		reference = new SeqDB(refNames, refSeqs, refAbunds, false);
		if (!reference->isNucleo()) {
		    return results;
		}
	}

	vector<ChimeHit2> Hits;
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

		string chimeraStatus = "N";
		if (Hit.Accept()) {  chimeraStatus = "Y"; }

		dataUchime->addOutput(Hit.Score, Hit.QLabel, Hit.ALabel, Hit.BLabel,
                                chimeraStatus, Hit.PctIdQM, Hit.PctIdQA,
                                Hit.PctIdQB, Hit.PctIdAB, Hit.PctIdQT,
                                Hit.Div, Hit.CS_LY, Hit.CS_LN, Hit.CS_LA,
                                Hit.CS_RY, Hit.CS_RN, Hit.CS_RA);
	}

	resultsNames.push_back("uchimeout");
	results.push_back(dataUchime->getUchimeOutput());
	results.attr("names") = resultsNames;

	return results;
}
/******************************************************************************/
