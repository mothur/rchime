/*
 *  chimera_uchime.cpp
 *
 *
 *  Created by Sarah Westcott on 2/20/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "chimera_uchime.h"
#include "chime.h"
#include "mothur.h"
#include "uchime_main.h"

/******************************************************************************/
ChimeraUchime::ChimeraUchime(bool derep, int proc, bool si, bool hg, Rcpp::List options) :
    Chimera(derep, proc, si, hg, options) {
}
/******************************************************************************/
void processUchime(chimeraData* params) {

    UchimeMain* uchime = new UchimeMain();

    if (params->groups.size() != 0) {

        SEXP* bar = nullptr;

        // main thread will update progress
        if (!params->silent) {
            // create progress bar
            int numGroups = params->groups.size();
            bar = new SEXP(PROTECT(cli_progress_bar(numGroups, NULL)));
            cli_progress_set_name(*bar, "Checking for chimeric sequences");
        }

        // for each group in dataset
        for (int k = 0; k < params->groups.size(); k++) {

            RcppThread::checkUserInterrupt();

            // main thread updates progress
            if (!params->silent) {
                cli_progress_set(*bar, k+1); 
            }

            // call uchime source code
            params->uchimeResults.push_back(uchime->runUchime(params->names[k], params->seqs[k],
                                                params->refNames, params->refSeqs,
                                                params->abunds[k], params->chimeras[k], &params->options));
        }

        // terminate progress bar
        if (!params->silent) {
            cli_progress_done(*bar);
            UNPROTECT(1);
            delete bar;
        }
    }else {
        params->uchimeResults.push_back(uchime->runUchime(params->names[0], params->seqs[0],
                                        params->refNames, params->refSeqs,
                                        params->abunds[0], params->chimeras[0], &params->options));
    }

    delete uchime;

}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset) {

    vector<ChimeHit2> uchimeOutputs;
    vector<string> resultsNames;

    if (hasGroupData) {

        Rcpp::Function getGroups = dataset["get_groups"];
        groups = Rcpp::as<vector<string>>(getGroups());

        // run uchime on each sample individually
        uchimeOutputs = createProcesses(dataset);

    } else {
        chimeraData* dataBundle = new chimeraData(dataset, groups, *opts);
        // run as one sample
        processUchime(dataBundle);

        uchimeOutputs = dataBundle->uchimeResults[0];
        vector<string> chimeras = toVector(dataBundle->chimeras[0]);

        delete dataBundle;
        removeChimerasFromDataset(dataset, chimeras, "chimera_uchime");
    }

    Rcpp::List results = createUchimeResults(uchimeOutputs);

    return results;
}
/******************************************************************************/
vector<ChimeHit2> ChimeraUchime::createProcesses(Rcpp::Environment& dataset) {
    vector<ChimeHit2> results;

    Utils util;
    vector<pieceOfWork> indexes = util.divideWork(groups.size(), processors);

    vector<RcppThread::Thread*> workerThreads;
    vector<chimeraData*> data;

    //Launch worker threads
    for (int i = 0; i < processors-1; i++) {
        // fill seqs with this processors groups
        vector<string> theseGroups;
        for (int j = indexes[i+1].start; j < indexes[i+1].end; j++) {
            theseGroups.push_back(groups[j]);
        }
        chimeraData* dataBundle = new chimeraData(dataset, theseGroups, *opts);
        data.push_back(dataBundle);
        workerThreads.push_back(new RcppThread::Thread(processUchime,
                                                       dataBundle));
    }

    // fill seqs with this processors groups
    vector<string> theseGroups;
    for (int j = indexes[0].start; j < indexes[0].end; j++) {
        theseGroups.push_back(groups[j]);
    }

    // only the main thread reports progress if silent is false
    chimeraData* dataBundle = new chimeraData(dataset, theseGroups, *opts, silent);
    processUchime(dataBundle);

    map<string, uchimeAbunds > seqAbunds = combineResults(dataBundle, data,
                                                         workerThreads, results);

    if (dereplicate) {
        // set dataset abundances
        Rcpp::Function setAbundances = dataset["set_abundances"];
        Rcpp::List abundances;
        vector<string> names;
        for (auto it = seqAbunds.begin(); it != seqAbunds.end(); it++) {
            names.push_back(it->first);
            abundances.push_back(it->second.abundances);
        }

        setAbundances(names, abundances, Rcpp::Named("reason",
                                                 "chimera_uchime"));
    }else {
        // seqAbunds just contains seqs to be removed
        vector<string> seqsToRemove;
        for (auto it = seqAbunds.begin(); it != seqAbunds.end(); it++) {
            seqsToRemove.push_back(it->first);
        }
        removeChimerasFromDataset(dataset, seqsToRemove, "chimera_uchime");
    }

    return results;
}
/******************************************************************************/
Rcpp::List ChimeraUchime::removeChimeras(Rcpp::Environment& dataset,
                                         Rcpp::Environment& reference) {

    chimeraData* dataBundle = new chimeraData(dataset, reference, *opts);

    // run as one sample
    processUchime(dataBundle);

    vector<ChimeHit2> hits = (dataBundle->uchimeResults[0]);
    vector<string> chimeras = toVector(dataBundle->chimeras[0]);
    delete dataBundle;

    removeChimerasFromDataset(dataset, chimeras, "chimera_uchime");

    Rcpp::List results = createUchimeResults(hits);

    return results;
}
/******************************************************************************/
// this is called when dereplicate = TRUE, meaning sequence are only removed
// from the samples they are flagged as chimeric in
map<string, uchimeAbunds > ChimeraUchime::combineResults(chimeraData*& dataBundle,
                                         vector<chimeraData*>& data,
                                         vector<RcppThread::Thread*>& workerThreads,
                                         vector<ChimeHit2>& results) {

    map<string, uchimeAbunds > seqAbunds;
    int groupIndex = 0;
    vector<vector<ChimeHit2> > uchimeOut = dataBundle->uchimeResults;

    // dereplicate = TRUE, guilty until proven innocent (all samples must find it to be chimeric)
    // dereplicate = FALSE, innocent until proven guilty (any sample can find it to be chimeric)

    // collect main threads results
    for (int i = 0; i < dataBundle->groups.size(); i++) {

        for (int j = 0; j < dataBundle->names[i].size(); j++) {

            string seqName = dataBundle->names[i][j];
            auto it = seqAbunds.find(seqName);

            // create new abunds vector and fill this groups abunds
            if (it == seqAbunds.end()) {
                
                vector<int> abunds(groups.size(), 0);
                uchimeAbunds abundFlag(abunds, dereplicate);

                if (dereplicate) {
                    
                    // if not chimeric in this sample, then add abundance
                    if (dataBundle->chimeras[i].count(seqName) == 0){
                        abundFlag.abundances[groupIndex] = dataBundle->abunds[i][j];
                        abundFlag.chimeric = false;
                    }
                    seqAbunds[seqName] = abundFlag;
                }else {
                    // if chimeric in this sample, then remove from dataset
                    if (dataBundle->chimeras[i].count(seqName) != 0){
                       abundFlag.chimeric = true;
                       seqAbunds[seqName] = abundFlag;
                    }
                } 
            }else{
                if (dereplicate) {
                    // if not chimeric in this sample, then add abundance
                    if (dataBundle->chimeras[i].count(seqName) == 0) {
                        it->second.abundances[groupIndex] = dataBundle->abunds[i][j];
                        it->second.chimeric = false;
                    }
                }
            }
        }
        groupIndex++;
    }

    delete dataBundle;

    // collect child thread results and join
    for (int i = 0; i < processors-1; i++) {

        workerThreads[i]->join();

        for (int j = 0; j < data[i]->groups.size(); j++) {

            uchimeOut.push_back(data[i]->uchimeResults[j]);

            for (int k = 0; k < data[i]->names[j].size(); k++) {

                string seqName = data[i]->names[j][k];
                auto it = seqAbunds.find(seqName);

                // create new abunds vector and fill this groups abunds
                if (it == seqAbunds.end()) {
                
                    vector<int> abunds(groups.size(), 0);
                    uchimeAbunds abundFlag(abunds, dereplicate);

                    if (dereplicate) {
                        
                        // if not chimeric in this sample, then add abundance
                        if (data[i]->chimeras[j].count(seqName) == 0){
                            abundFlag.abundances[groupIndex] = data[i]->abunds[j][k];
                            abundFlag.chimeric = false;
                        }
                        seqAbunds[seqName] = abundFlag;
                    }else {
                        // if chimeric in this sample, then remove from dataset
                        if (data[i]->chimeras[j].count(seqName) != 0){
                            abundFlag.chimeric = true;
                            seqAbunds[seqName] = abundFlag;
                        }
                    }
                }else{
                    if (dereplicate) {
                        // if not chimeric in this sample, then add abundance
                        if (data[i]->chimeras[j].count(seqName) == 0) {
                            it->second.abundances[groupIndex] = data[i]->abunds[j][k];
                            it->second.chimeric = false;
                        }
                    }
                }
            }
            groupIndex++;
        }

        delete data[i];
        delete workerThreads[i];
    }

    // merge uchime results using flagged chimeras
    // the result will have 1 entry in hits for each sequence in dataset

    // if all samples find the query to be chimeric, choose first chimeric result
    // else pick first non-chimeric result.

    double minh = opts->getMinh();
	double mindiv = opts->getMindiv();
	int mindiffs = opts->getMindiffs();
    set<string> resolved;
    for (int i = 0; i < uchimeOut.size(); i++) {

        for (int j = 0; j < uchimeOut[i].size(); j++) {

            string seqName = uchimeOut[i][j].QLabel;
            auto it = resolved.find(seqName);

            // if we haven't recorded this sequence
            if (it == resolved.end()) {

                // dereplicate == true all seqs are in seqAbunds
                // dereplicate == false only chimeric seqs are in seqAbunds
                if (dereplicate) {

                    // sequence was found to be chimeric in all samples
                    if (seqAbunds[seqName].chimeric) {
                        results.push_back(uchimeOut[i][j]);
                        resolved.insert(seqName);
                    }else {
                        // pick first non-chimeric results
                        if (!uchimeOut[i][j].Accept(minh, mindiv, mindiffs)) {
                            results.push_back(uchimeOut[i][j]);
                            resolved.insert(seqName);
                        }
                    }
                }else {
                    // if chimeric, then save
                    if (uchimeOut[i][j].Accept(minh, mindiv, mindiffs)) {
                        results.push_back(uchimeOut[i][j]);
                        resolved.insert(seqName);
                    }else {
                        auto itChime = seqAbunds.find(seqName);

                        // non chimera
                        if (itChime == seqAbunds.end()) {
                            results.push_back(uchimeOut[i][j]);
                            resolved.insert(seqName);
                        }
                    }  
                }
            }
        }
    }

    return seqAbunds;
}
/******************************************************************************/
// uchimealns
// Each chimeric sequence creates the following:
//
// Query   (  179 nt) F21Fcsw_11639
// ParentA (  179 nt) F11Fcsw_6529
// ParentB (  181 nt) F21Fcsw_12128

// A     1 AAGgAAGAtTAATACaagATGgCaTCatgAGtccgCATgTtcAcatGATTAAAG--gTaTtcCGGTagacGATGGGGATG 78
// Q     1 AAGTAAGACTAATACCCAATGACGTCTCTAGAAGACATCTGAAAGAGATTAAAG--ATTTATCGGTGATGGATGGGGATG 78
// B     1 AAGgAAGAtTAATcCaggATGggaTCatgAGttcACATgTccgcatGATTAAAGgtATTTtcCGGTagacGATGGGGATG 80
// Diffs      N    N    A N?N   N N  NNN  N?NB   N ?NaNNN          B B NN    NNNN
// Votes      0    0    + 000   0 0  000  000+   0 00!000            + 00    0000
// Model   AAAAAAAAAAAAAAAAAAAAAAxBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

// A    79 CGTtccATTAGaTaGTaGGCGGGGTAACGGCCCACCtAGtCttCGATggaTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 158
// Q    79 CGTCTGATTAGCTTGTTGGCGGGGTAACGGCCCACCAAGGCAACGATCAGTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 158
// B    81 CGTtccATTAGaTaGTaGGCGGGGTAACGGCCCACCtAGtCAACGATggaTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 160
// Diffs      NNN     N N  N                   N  N BB    NNN
// Votes      000     0 0  0                   0  0 ++    000
// Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

// A   159 TGGAACTGAGACACGGTCCAA 179
// Q   159 TGGAACTGAGACACGGTCCAA 179
// B   161 TGGAACTGAGACACGGTCCAA 181
// Diffs
// Votes
// Model   BBBBBBBBBBBBBBBBBBBBB

// Ids.  QA 76.6%, QB 77.7%, AB 93.7%, QModel 78.9%, Div. +1.5%
// Diffs Left 7: N 0, A 6, Y 1 (14.3%); Right 35: N 1, A 30, Y 4 (11.4%), Score 0.0047
/******************************************************************************/
Rcpp::List ChimeraUchime::createUchimeResults(vector<ChimeHit2> hits) {

    Rcpp::List results;
    vector<string> resultsNames;

    // Q, A, B, Y
    vector<string> queries, Aparents, Bparents, chimericStatus;

    // Score(Higher score indicates chimeric status), Divs(IdQM - IdQT)
    // IdQM, IdQA, IdQB, IdQT, IdAB - similarity scores
    vector<double> scores, IdQMs, IdQAs, IdQBs, IdABs, IdQTs, Divs;

    // LY, LN, LA, RN, RY, RA - segment votes
    vector<int> LYs, LNs, LAs, RNs, RYs, RAs;

    Rcpp::List uchimeAlns;
    vector<string> chimeras;

    double minh = opts->getMinh();
	double mindiv = opts->getMindiv();
	int mindiffs = opts->getMindiffs();

    // create a uchimeout
    for (int i = 0; i < hits.size(); i++) {

        string status = "N";
        if (hits[i].Accept(minh, mindiv, mindiffs)) {
            status = "Y";
            chimeras.push_back(hits[i].QLabel);

            //create uchimealns record
            if (opts->getChimealns()) {

                string query = hits[i].Q3;
                string parentA = hits[i].A3;
                string parentB = hits[i].B3;
                int length = query.length();

                vector<string> records;
                records.push_back("Query ( " + toString(getNumBases(query)) + "nt ) " + hits[i].QLabel); //+ "\n";
                records.push_back("ParentA ( " + toString(getNumBases(parentA)) + "nt ) " + hits[i].ALabel); // + "\n";
                records.push_back("ParentB ( " + toString(getNumBases(parentB)) + "nt ) " + hits[i].BLabel); // + "\n\n";

                int firstBase, lastBase, qPos, aPos, bPos;
                firstBase = MOTHURMAX;
                aPos = 0; bPos = 0; qPos = 0;
                // Strip terminal gaps in query
                for (int j = 0; j < length; j++) {
                    // if not a gap
                    if (!isGap(query[j])) {
                        if (firstBase == MOTHURMAX) {
                            firstBase = j;
                        }
                        lastBase = j;
                    }
                }

                // set parents positions
                for (int j = 0; j < firstBase; j++) {
                    // if not a gap
                    if (!isGap(parentA[j])) { aPos++; }
                    if (!isGap(parentB[j])) { bPos++; }
                }

                // create Q, A, B, Diffs, Votes, and Model rows
                string A, Q, B, Diffs, Votes, Model;
                A = "A " + toString(aPos+1) + " ";
                Q = "Q " + toString(qPos+1) + " ";
                B = "B " + toString(bPos+1) + " ";
                Diffs = "Diffs   ";
                Votes = "Votes   ";
                Model = "Model   ";
                for (int j = firstBase; j <= lastBase; j++) {

                    char q = query[j];
                    char a = parentA[j];
                    char b = parentB[j];

                    // Q, A, B rows
                    Q += q;
                    if (a != q) { a = tolower(a); }
                    if (b != q) { b = tolower(b); }
                    A += a;
                    B += b;
                    if ((a != '.') && (a != '-')) { aPos++; }
                    if ((b != '.') && (b != '-')) { bPos++; }

                    // Diffs row
                    char c = ' ';
                    if (isGap(q) || isGap(a) || isGap(b)) {
                        c = ' ';
                    }else if (j < hits[i].ColXLo) {
				        if (q == a && q == b)      { c = ' '; }
                        else if (q == a && q != b) { c = 'A'; }
                        else if (q == b && q != a) { c = 'b'; }
                        else if (a == b && q != a) { c = 'N'; }
                        else                       { c = '?'; }
				    }else if (j > hits[i].ColXHi) {
				        if (q == a && q == b)      { c = ' '; }
                        else if (q == b && q != a) { c = 'B'; }
                        else if (q == a && q != b) { c = 'a'; }
                        else if (a == b && q != a) { c = 'N'; }
                        else                       { c = '?'; }
				    }
                    Diffs += c;

                    // Votes row
                    bool PrevGap = false;
                    if (j > 0) {
                        PrevGap = (isgap(query[j-1]) || isgap(parentA[j-1]) || isgap(parentB[j-1]));
                    }
                    bool NextGap = false;
                    if (j+1 < length) {
                        NextGap = (isgap(query[j+1]) || isgap(parentA[j+1]) || isgap(parentB[j+1]));
                    }

                    c = ' ';
			        if (isgap(q) || isgap(a) || isgap(b) || PrevGap || NextGap) { }
                    else if (j < hits[i].ColXLo) {
				        if (q == a && q == b)        { c = ' '; }
                        else if (q == a && q != b)   { c = '+'; }
                        else if (q == b && q != a)   { c = '!'; }
                        else                         { c = '0'; }
			        }else if (j > hits[i].ColXHi) {
				        if (q == a && q == b)        { c = ' '; }
				        else if (q == b && q != a)   { c = '+'; }
				        else if (q == a && q != b)   { c = '!'; }
				        else                         { c = '0'; }
			        }
                    Votes += c;

                    // Model row
                    if (j < hits[i].ColXLo)   { Model += "A"; }
                    else if (j >= hits[i].ColXLo && j <= hits[i].ColXHi) { Model += "x"; }
                    else                  { Model += "B"; }


                }

                records.push_back(A); //+ "\n";
                records.push_back(Q); // + "\n";
                records.push_back(B); // + "\n";
                records.push_back(Diffs); // + "\n";
                records.push_back(Votes); // + "\n";
                records.push_back(Model); // + "\n\n";

                // Ids row
	            double PctIdBestP = max(hits[i].PctIdQA, hits[i].PctIdQB);
	            double Div = (hits[i].PctIdQM - PctIdBestP)*100.0/PctIdBestP;

	            unsigned LTot = hits[i].CS_LY + hits[i].CS_LN + hits[i].CS_LA;
	            unsigned RTot = hits[i].CS_RY + hits[i].CS_RN + hits[i].CS_RA;

                string idsRow = "Ids.  QA " + toString(hits[i].PctIdQA);
                idsRow += ", QB " + toString(hits[i].PctIdQB);
                idsRow += ", AB " + toString(hits[i].PctIdAB);
                idsRow += ", QModel " + toString(hits[i].PctIdQM);
                idsRow += ", Div. " + toString(Div); // + "\n";
                records.push_back(idsRow);

                double PctL = Pct(hits[i].CS_LY, LTot);
	            double PctR = Pct(hits[i].CS_RY, RTot);

                // Diffs row
                string diffsRow = "Diffs Left " + toString(LTot) + ": N " +  toString(hits[i].CS_LN);
                diffsRow += ", A " + toString(hits[i].CS_LA) + ", Y " + toString(hits[i].CS_LY);
                diffsRow += "(" + toString(PctL) + ");";
                diffsRow += " Right " + toString(RTot) + ": N " +  toString(hits[i].CS_RN);
                diffsRow += ", A " + toString(hits[i].CS_RA) + ", Y " + toString(hits[i].CS_RY);
                diffsRow += "(" + toString(PctR) + "), Score " + toString(hits[i].Score); // + "\n";
                records.push_back(diffsRow);

                uchimeAlns.push_back(records);
            }
        }

        chimericStatus.push_back(status);
        
        // uchimeout
        if (hits[i].Div <= 0.0) {
            scores.push_back(0.0000);
            queries.push_back(hits[i].QLabel+"/ab="+toString(hits[i].AbQ)+"/");
            Aparents.push_back("*");
            Bparents.push_back("*");
            IdQMs.push_back(-1.0);
            IdQAs.push_back(-1.0);
            IdQBs.push_back(-1.0);
            IdABs.push_back(-1.0);
            IdQTs.push_back(-1.0);
            Divs.push_back(-1.0);
            LYs.push_back(-1);
            LNs.push_back(-1);
            LAs.push_back(-1);
            RNs.push_back(-1);
            RYs.push_back(-1);
            RAs.push_back(-1);
        }else {
            scores.push_back(hits[i].Score);
            queries.push_back(hits[i].QLabel+"/ab="+toString(hits[i].AbQ)+"/");
            Aparents.push_back(hits[i].ALabel+"/ab="+toString(hits[i].AbA)+"/");
            Bparents.push_back(hits[i].BLabel+"/ab="+toString(hits[i].AbB)+"/");
            IdQMs.push_back(hits[i].PctIdQM);
            IdQAs.push_back(hits[i].PctIdQA);
            IdQBs.push_back(hits[i].PctIdQB);
            IdABs.push_back(hits[i].PctIdAB);
            IdQTs.push_back(hits[i].PctIdQT);
            Divs.push_back(hits[i].Div);
            LYs.push_back(hits[i].CS_LY);
            LNs.push_back(hits[i].CS_LN);
            LAs.push_back(hits[i].CS_LA);
            RNs.push_back(hits[i].CS_RN);
            RYs.push_back(hits[i].CS_RY);
            RAs.push_back(hits[i].CS_RA);
        }
    }

    Rcpp::DataFrame uchimeOut = DataFrame::create(
         Named("Score") = scores,
         _("Query") = queries, _("A") = Aparents, _("B") = Bparents,
         _("IdQM") = IdQMs, _("IdQA") = IdQAs, _("IdQB") = IdQBs,
         _("IdAB") = IdABs, _("IdQT") = IdQTs,
         _("LY") = LYs, _("LN") = LNs, _("LA") = LAs,
         _("RY") = RYs, _("RN") = RNs, _("RA") = RAs,
         _("Div") = Divs, _("Y") = chimericStatus);

    results.push_back(uchimeOut);
    resultsNames.push_back("uchimeout");

    results.push_back(chimeras);
    resultsNames.push_back("accnos");

    if (opts->getChimealns()) {
        resultsNames.push_back("uchimealns");
        results.push_back(uchimeAlns);
    }
    results.attr("names") = resultsNames;

    return results;
}
/******************************************************************************/
double ChimeraUchime::Pct(double x, double y) {
    if (y == 0.0f) {
        return 0.0f;
    }
    return (x*100.0f)/y;
}
/******************************************************************************/
