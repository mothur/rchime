
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
        opts->opt_minwordmatches = minwordmatches_defaults[opts->opt_wordlength];
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
        sortedVector[i].abund = (abunds)[i];
    }
    sort(sortedVector.begin(), sortedVector.end(), compareAbundance);

    vector<unsigned> order((abunds).size(), 0);
    for (int i = 0; i < (abunds).size(); i++) {
        order[i] = sortedVector[i].index;
        (abunds)[i] = sortedVector[i].abund;
    }

    applyOrder(sequenceNames, order);
    applyOrder(sequences, order);
    applyOrder(abunds, order);
}
/******************************************************************************/
void processDenovoVsearch(chimeraData* params) {

    std::vector<std::string> pRefNames;
    std::vector<std::string> pRefSequences;

    // denovo by sample
    for (int k = params->start; k < params->stop; k++) {

        RcppThread::checkUserInterrupt();

        Vsearch_Main vsearch;

        // call vsearch source code
        params->results[k] = vsearch.vmain(params->pNames->at(k),
                                            params->pSequences->at(k),
                                            pRefNames, pRefSequences,
                                            params->pAbundances->at(k), 0, 0);

        params->chimeras.push_back(vsearch.chimeraNames);
    }
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

    // denovo single sample
    if (!bySample) {
        // blank reference
        std::vector<std::string> refNames;
        std::vector<std::string> refSeqs;

        Vsearch_Main vsearch;

        sortDescending(pNames->at(0),
                       pSequences->at(0),
                       pAbundances->at(0));

        vector<ChimeHit2> vsearchOutput = vsearch.vmain(pNames->at(0),
                                                         pSequences->at(0),
                                                         refNames, refSeqs,
                                                         pAbundances->at(0),
                                                         0, 0);

        return createVsearchResults(vsearchOutput);
    }

    // denovo - run vsearch on each sample individually
    return createDenovoProcesses();
}
/******************************************************************************/
Rcpp::List ChimeraVsearch::createDenovoProcesses() {

    vector<RcppThread::Thread*> workerThreads;
    vector<chimeraData*> data;

    // divides samples between processors
    vector<int> starts, ends;
    vector<pieceOfWork> indexes = divideWork(pNames->size(), processors);

    //Launch worker threads
    for (int i = 0; i < processors-1; i++) {

        chimeraData* dataBundle = new chimeraData(pNames, pSequences, pAbundances,
                                                  indexes[i+1].start,
                                                  indexes[i+1].end);
        data.push_back(dataBundle);
        workerThreads.push_back(new RcppThread::Thread(processDenovoVsearch,
                                                       dataBundle));
    }

    chimeraData* dataBundle = new chimeraData(pNames, pSequences, pAbundances,
                                              indexes[0].start,
                                              indexes[0].end);
    processDenovoVsearch(dataBundle);

    vector<ChimeHit2> vsearchResults;
    map<string, vsearchAbunds> seqAbunds = combineResults(dataBundle,
                                                           data,
                                                           workerThreads,
                                                           vsearchResults);

    Rcpp::List results = createVsearchResults(vsearchResults);

    if (dereplicate) {
        // add seqAbunds results to results so the sequences can be removed
        // by sample in rchime function
        Rcpp::List newAbundances(seqAbunds.size());

        int i = 0;
        for (const auto& pair : seqAbunds) {
            newAbundances[i] = pair.second.abundances;
            i++;
        }

        newAbundances.names() = getKeys(seqAbunds);

        vector<string> resultsNames = results.names();
        resultsNames.push_back("set_abundance_values");

        // add map to results
        results.push_back(newAbundances);
        results.attr("names") = resultsNames;
    }

    return results;
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
// this is called when dereplicate = TRUE, meaning sequence are only removed
// from the samples they are flagged as chimeric in
map<string, vsearchAbunds > ChimeraVsearch::combineResults(chimeraData*& dataBundle,
                                         vector<chimeraData*>& data,
                                         vector<RcppThread::Thread*>& workerThreads,
                                         vector<ChimeHit2>& results) {

    map<string, vsearchAbunds > seqAbunds;
    int groupIndex = dataBundle->start;
    vector<vector<ChimeHit2> > vsearchOut = dataBundle->results;

    // dereplicate = TRUE, guilty until proven innocent (all samples must find it to be chimeric)
    // dereplicate = FALSE, innocent until proven guilty (any sample can find it to be chimeric)

    // collect main threads results
    for (int i = dataBundle->start; i < dataBundle->stop; i++) {

        for (int j = 0; j < dataBundle->pNames->at(i).size(); j++) {

            string seqName = dataBundle->pNames->at(i)[j];
            auto it = seqAbunds.find(seqName);

            // create new abunds vector and fill this groups abunds
            if (it == seqAbunds.end()) {

                vector<float> abunds((dataBundle->pNames->size()), 0);
                vsearchAbunds abundFlag(abunds, dereplicate);

                if (dereplicate) {

                    // if not chimeric in this sample, then add abundance
                    if (dataBundle->chimeras[i].count(seqName) == 0){
                        abundFlag.abundances[groupIndex] = dataBundle->pAbundances->at(i)[j];
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
                        it->second.abundances[groupIndex] = dataBundle->pAbundances->at(i)[j];
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

        for (int j = data[i]->start; j < data[i]->stop; j++) {

            vsearchOut.push_back(data[i]->results[j]);

            for (int k = 0; k < data[i]->pNames->at(j).size(); k++) {

                string seqName = data[i]->pNames->at(j)[k];
                auto it = seqAbunds.find(seqName);

                // create new abunds vector and fill this groups abunds
                if (it == seqAbunds.end()) {

                    vector<float> abunds((dataBundle->pNames->size()), 0);
                    vsearchAbunds abundFlag(abunds, dereplicate);

                    if (dereplicate) {

                        // if not chimeric in this sample, then add abundance
                        if (data[i]->chimeras[j].count(seqName) == 0){
                            abundFlag.abundances[groupIndex] = data[i]->pAbundances->at(j)[k];
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
                            it->second.abundances[groupIndex] = data[i]->pAbundances->at(j)[k];
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

    // merge vsearch results using flagged chimeras
    // the result will have 1 entry in hits for each sequence in dataset

    // if all samples find the query to be chimeric, choose first chimeric result
    // else pick first non-chimeric result.
    set<string> resolved;
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
                        if (vsearchOut[i][j].status == "N") {
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

    return seqAbunds;
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

    Rcpp::DataFrame uchimeOut = Rcpp::DataFrame::create(
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

    results.push_back(uchimeOut);
    resultsNames.push_back("uchimeout");

    results.push_back(chimeras);
    resultsNames.push_back("accnos");

    results.attr("names") = resultsNames;

    return results;
}
/******************************************************************************/
