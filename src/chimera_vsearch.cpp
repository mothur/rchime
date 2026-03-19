
#include "chimera_vsearch.h"
#include "vsearch_main.h"

Vsearch_Options* Vsearch_Options::_uniqueInstance = 0;

/******************************************************************************/
ChimeraVsearch::ChimeraVsearch(std::vector<std::vector<string>>& sequenceNames,
                               std::vector<std::vector<string>>& sequences,
                               std::vector<std::vector<float>>& abundances,
                               Rcpp::List options)  {

    pRefNames = nullptr;
    pRefSequences = nullptr;
    pNames = &sequenceNames;
    pSequences = &sequences;
    pAbundances = &abundances;
    opts = Vsearch_Options::getInstance();
    opts->reset();

    bySample = false;
    if ((*pNames).size() > 1) {
        bySample = true;
    }

    Rcpp::CharacterVector optionsNames = options.names();

    if (contains("abskew", optionsNames))   {
        opts->opt_abskew = Rcpp::as<double>(options["abskew"]); }
    if (contains("minh", optionsNames))   {
        opts->opt_minh = Rcpp::as<double>(options["minh"]);   }
    if (contains("mindiv", optionsNames))   {
        opts->opt_mindiv = Rcpp::as<double>(options["mindiv"]);   }
    if (contains("xn", optionsNames))   {
        opts->opt_xn = Rcpp::as<double>(options["xn"]);   }
    if (contains("dn", optionsNames))   {
        opts->opt_dn = Rcpp::as<float>(options["dn"]);   }

    opts->opt_threads = Rcpp::as<int>(options["processors"]);
    processors = opts->opt_threads;
    dereplicate = Rcpp::as<bool>(options["dereplicate"]);

    // set once - so moved out of vsearch code
    opts->opt_gap_open_query_left -= opts->opt_gap_extension_query_left;
    opts->opt_gap_open_target_left -= opts->opt_gap_extension_target_left;
    opts->opt_gap_open_query_interior -= opts->opt_gap_extension_query_interior;
    opts->opt_gap_open_target_interior -= opts->opt_gap_extension_target_interior;
    opts->opt_gap_open_query_right -= opts->opt_gap_extension_query_right;
    opts->opt_gap_open_target_right -= opts->opt_gap_extension_target_right;

    if (opts->opt_minwordmatches < 0) {
        opts->opt_minwordmatches = minwordmatches_defaults[8];
    }

}
/******************************************************************************/
bool ChimeraVsearch::contains(std::string s, Rcpp::CharacterVector& nv) {
    for (int i=0; i<nv.size(); i++) {
        if (std::string(nv[i]) == s) {
            return true;
        }
    }
    return false;
}
/******************************************************************************/
void ChimeraVsearch::sortDescending(std::vector<std::string>& sequenceNames,
                                    std::vector<std::string>& sequences,
                                    std::vector<float>& abunds) {

    vector<orderFloatAbundance> sortedVector((abunds).size());

    for (int i = 0; i < (abunds).size(); i++) {
        sortedVector[i].index = i;
        sortedVector[i].abund = abunds[i];
        sortedVector[i].name = sequenceNames[i];
    }

    sort(sortedVector.begin(), sortedVector.end(), compareAbundance);

    vector<unsigned> order((abunds).size(), 0);
    for (int i = 0; i < (abunds).size(); i++) {
        order[i] = sortedVector[i].index;
        abunds[i] = sortedVector[i].abund;
        sequenceNames[i] = sortedVector[i].name;
    }

    applyOrder(sequences, order);
}
/******************************************************************************/
void processReferenceVsearch(chimeraData* params) {

    Vsearch_Main vsearch;

    // call vsearch source code
    params->results[0] = vsearch.vmain(params->pNames->at(0),
                                       params->pSequences->at(0),
                                       *params->pRefNames, *params->pRefSequences,
                                       params->pAbundances->at(0),
                                       params->start, params->stop);

}
/******************************************************************************/
Rcpp::List ChimeraVsearch::removeChimeras() {

    // blank reference
    std::vector<std::string> refNames;
    std::vector<std::string> refSeqs;
    Rcpp::List denovoResults;

    // denovo single sample
    if (!bySample) {

        Vsearch_Main vsearch;

        sortDescending(pNames->at(0),
                       pSequences->at(0),
                       pAbundances->at(0));

        vector<ChimeHit2> results = vsearch.vmain(pNames->at(0),
                                pSequences->at(0),
                                refNames, refSeqs,
                                pAbundances->at(0), 0, 0);

        denovoResults = createVsearchResults(results);
    }else{

        // outputs
        vector<vector<ChimeHit2>> vsearchOut(pNames->size());
        map<string, vsearchAbunds > seqAbunds;

        // denovo by sample
        for (int k = 0; k < pNames->size(); k++) {

            RcppThread::checkUserInterrupt();

            sortDescending(pNames->at(k),
                           pSequences->at(k),
                           pAbundances->at(k));

            Vsearch_Main vsearch;

            // call vsearch source code
            vsearchOut[k] = vsearch.vmain(pNames->at(k),
                                        pSequences->at(k),
                                        refNames, refSeqs,
                                        pAbundances->at(k), 0, 0);

            set<string> chimeras = vsearch.chimeraNames;

            // combining results
            for (int j = 0; j < pNames->at(k).size(); j++) {

                string seqName = pNames->at(k)[j];
                auto it = seqAbunds.find(seqName);

                // create new abunds vector and fill this groups abunds
                if (it == seqAbunds.end()) {

                    // abunds(numSamples, 0)
                    vector<float> abunds((pNames->size()), 0);
                    vsearchAbunds abundFlag(abunds, dereplicate);

                    if (dereplicate) {

                        // if not chimeric in this sample, then add abundance
                        if (chimeras.count(seqName) == 0){
                            abundFlag.abundances[k] = pAbundances->at(k)[j];
                            abundFlag.chimeric = false;
                        }
                        seqAbunds[seqName] = abundFlag;
                    }else {
                        // if chimeric in this sample, then remove from dataset
                        if (chimeras.count(seqName) != 0){
                            abundFlag.chimeric = true;
                            seqAbunds[seqName] = abundFlag;
                        }
                    }
                }else{
                    if (dereplicate) {
                        // if not chimeric in this sample, then add abundance
                        if (chimeras.count(seqName) == 0) {
                            it->second.abundances[k] = pAbundances->at(k)[j];
                            it->second.chimeric = false;
                        }
                    }
                }
            }

            pNames->at(k).clear();
            pSequences->at(k).clear();
            pAbundances->at(k).clear();
        }

        // merge vsearch results using flagged chimeras
        // the result will have 1 entry in hits for each sequence in dataset

        // if all samples find the query to be chimeric, choose first chimeric result
        // else pick first non-chimeric result.
        set<string> resolved;
        vector<ChimeHit2> results;
        for (int i = 0; i < vsearchOut.size(); i++) {

            for (int j = 0; j < vsearchOut[i].size(); j++) {

                string seqName = vsearchOut[i][j].QLabel;
                auto it = resolved.find(seqName);

                // if we haven't recorded this sequence
                if (it == resolved.end()) {

                    // dereplicate == true all seqs are in seqAbunds
                    // dereplicate == false only chimeric seqs are in seqAbunds
                    if (dereplicate) {

                        // sequence was found to be chimeric in all samples
                        if (seqAbunds[seqName].chimeric) {
                            results.push_back(vsearchOut[i][j]);
                            resolved.insert(seqName);
                        }else {
                            // pick first non-chimeric results
                            if (vsearchOut[i][j].status != "Y") {
                                results.push_back(vsearchOut[i][j]);
                                resolved.insert(seqName);
                            }
                        }
                    }else {
                        // if chimeric, then save
                        if (vsearchOut[i][j].status == "Y") {
                            results.push_back(vsearchOut[i][j]);
                            resolved.insert(seqName);
                        }else {
                            auto itChime = seqAbunds.find(seqName);

                            // non chimera
                            if (itChime == seqAbunds.end()) {
                                results.push_back(vsearchOut[i][j]);
                                resolved.insert(seqName);
                            }
                        }
                    }
                }
            }
        }

        denovoResults = createVsearchResults(results);

        if (dereplicate) {
            // add seqAbunds results to results so the sequences can be removed
            // by sample in rchime function
            Rcpp::List newAbundances(seqAbunds.size());

            int i = 0;
            for (auto it = seqAbunds.begin(); it != seqAbunds.end(); it++) {
                newAbundances[i] = it->second.abundances;
                i++;
            }

            vector<string> resultsNames = denovoResults.names();
            resultsNames.push_back("set_abundance_values");

            Rcpp::List abundance_values = Rcpp::List::create(
                Rcpp::Named("names") = getKeys(seqAbunds),
                Rcpp::_("abundances") = newAbundances);

            // add map to results
            denovoResults.push_back(abundance_values);
            denovoResults.attr("names") = resultsNames;
        }

    }

    return denovoResults;
}
/******************************************************************************/
Rcpp::List ChimeraVsearch::removeChimeras(std::vector<std::string>& refNames,
                                          std::vector<std::string>& refSequences) {

    dereplicate = false;
    pRefNames = &refNames;
    pRefSequences = &refSequences;

    vector<RcppThread::Thread*> workerThreads;
    vector<chimeraData*> data;

    vector<int> starts, ends;
    vector<pieceOfWork> indexes = divideWork(pNames->at(0).size(), processors);

    //Launch worker threads
    for (int i = 0; i < processors-1; i++) {

        chimeraData* dataBundle = new chimeraData(pNames, pSequences, pAbundances,
                                                  &refNames, &refSequences,
                                                  indexes[i+1].start,
                                                  indexes[i+1].end);
        data.push_back(dataBundle);
        workerThreads.push_back(new RcppThread::Thread(processReferenceVsearch,
                                                       dataBundle));
    }

    chimeraData* dataBundle = new chimeraData(pNames, pSequences, pAbundances,
                                              &refNames, &refSequences,
                                              indexes[0].start,
                                              indexes[0].end);
    processReferenceVsearch(dataBundle);

    vector<ChimeHit2> hits = dataBundle->results[0];

    for (int i = 0; i < workerThreads.size(); i++) {

        workerThreads[i]->join();

        hits.insert(hits.end(),
                    data[i]->results[0].begin(), data[i]->results[0].end());

        delete data[i];
        delete workerThreads[i];
    }

    Rcpp::List results = createVsearchResults(hits);
    return results;
}
/******************************************************************************/
Rcpp::List ChimeraVsearch::createVsearchResults(vector<ChimeHit2> hits) {

    Rcpp::List results;
    vector<string> resultsNames;

    // Default constructor values
    // QLabel = ""; ALabel = "*"; BLabel = "*"; CLabel = "*"; status = "N";
    // QM = 0.0; QA = 0.0; QB = 0.0; QC = 0.0; QT = 0.0;
    // LY = 0.0; LN = 0.0; LA = 0.0; RY = 0.0; RN = 0.0; RA = 0.0;
    // Div = 0.0; Score = 0.0; H = 0.0;

    // Q, A, B, Y
    vector<string> queries, Aparents, Bparents, Cparents, chimericStatus;
    vector<double> scores, QMs, QAs, QBs, QCs, QTs, Divs;

    // LY, LN, LA, RN, RY, RA - segment votes
    vector<int> LYs, LNs, LAs, RNs, RYs, RAs;
    vector<string> chimeras;

    // create a uchimeout
    for (int i = 0; i < hits.size(); i++) {

        if (hits[i].status == "Y") {
            chimeras.push_back(hits[i].QLabel);
        }

        chimericStatus.push_back(hits[i].status);
        queries.push_back(hits[i].QLabel);

        // uchimeout
        if (hits[i].Div <= 0.0) {
            scores.push_back(0.0000);
            Aparents.push_back("*");
            Bparents.push_back("*");
            Cparents.push_back("*");
            QMs.push_back(-1.0);
            QAs.push_back(-1.0);
            QBs.push_back(-1.0);
            QCs.push_back(-1.0);
            QTs.push_back(-1.0);
            Divs.push_back(-1.0);
            LYs.push_back(-1);
            LNs.push_back(-1);
            LAs.push_back(-1);
            RNs.push_back(-1);
            RYs.push_back(-1);
            RAs.push_back(-1);
        }else {
            scores.push_back(hits[i].Score);
            Aparents.push_back(hits[i].ALabel);
            Bparents.push_back(hits[i].BLabel);
            Cparents.push_back(hits[i].CLabel);
            QMs.push_back(hits[i].QM);
            QAs.push_back(hits[i].QA);
            QBs.push_back(hits[i].QB);
            QCs.push_back(hits[i].QC);
            QTs.push_back(hits[i].QT);
            LYs.push_back(hits[i].LY);
            LNs.push_back(hits[i].LN);
            LAs.push_back(hits[i].LA);
            RNs.push_back(hits[i].RN);
            RYs.push_back(hits[i].RY);
            RAs.push_back(hits[i].RA);
            Divs.push_back(hits[i].Div);
        }
    }

    Rcpp::DataFrame rchime_report = Rcpp::DataFrame::create(
        Rcpp::Named("Score") = scores,
        Rcpp::_("Query") = queries,
        Rcpp::_("ParentA") = Aparents,
        Rcpp::_("ParentB") = Bparents,
        Rcpp::_("ParentC") = Cparents,
        Rcpp::_("QM") = QMs, Rcpp::_("QA") = QAs, Rcpp::_("QB") = QBs,
        Rcpp::_("QC") = QCs, Rcpp::_("QT") = QTs,
        Rcpp::_("LY") = LYs, Rcpp::_("LN") = LNs, Rcpp::_("LA") = LAs,
        Rcpp::_("RY") = RYs, Rcpp::_("RN") = RNs, Rcpp::_("RA") = RAs,
        Rcpp::_("Div") = Divs, Rcpp::_("Chimeric_Status") = chimericStatus);

    results.push_back(rchime_report);
    resultsNames.push_back("rchime_report");

    results.push_back(chimeras);
    resultsNames.push_back("accnos");

    results.attr("names") = resultsNames;

    return results;
}
/******************************************************************************/
