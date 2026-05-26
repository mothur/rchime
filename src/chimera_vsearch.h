#ifndef CHIMERAVSEARCH_H
#define CHIMERAVSEARCH_H

#include "../inst/include/rchime.h"
#include "chime.h"
#include "options.h"

/******************************************************************************/
struct orderFloatAbundance {
    std::string name;
    float abund;
    unsigned index;

    orderFloatAbundance() { abund = 0; index = 0; name = ""; }
    orderFloatAbundance(std::string n, float a, unsigned i) : name(n), abund(a), index(i) {}
    ~orderFloatAbundance() {}
};
/******************************************************************************/
//sorts highest to lowest
static inline bool compareAbundance(orderFloatAbundance left, orderFloatAbundance right){
    if (left.abund > right.abund) {
        return true;
    }

    if (left.abund < right.abund) {
        return false;
    }

    // the abundance is equal, sort alphabetically
    if (std::strcmp(left.name.c_str(), right.name.c_str()) < 0) {
        return true;
    }
    return false;
}
/******************************************************************************/
struct vsearchAbunds {
    std::vector<float> abundances;
    bool chimeric;

    vsearchAbunds() { chimeric = true; }
    vsearchAbunds(std::vector<float> abunds, bool c) {
        abundances = abunds;
        chimeric = c;
    }
    ~vsearchAbunds() {}

};
/******************************************************************************/
struct chimeraData {

    // inputs
    std::vector<std::vector<std::string>>* pNames;
    std::vector<std::vector<std::string>>* pSequences;
    std::vector<std::vector<float>>* pAbundances;

    std::vector<std::string>* pRefNames;
    std::vector<std::string>* pRefSequences;

    int start;
    int stop;

    // outputs
    std::vector<std::vector<ChimeHit2>> results;

    // reference - each processor gets assigned part of names and sequences
    chimeraData(std::vector<std::vector<std::string>>* names,
                std::vector<std::vector<std::string>>* sequences,
                std::vector<std::vector<float>>* abundances,
                std::vector<std::string>* rnames,
                std::vector<std::string>* rsequences,
                int s, int e) {

        start = s;
        stop = e;

        pNames = names;
        pSequences = sequences;
        pAbundances = abundances;
        pRefNames = rnames;
        pRefSequences = rsequences;

        results.resize(1);
    }
    ~chimeraData() {}
};
/******************************************************************************/
class ChimeraVsearch {

public:

    // dereplicate, processors, silent
    ChimeraVsearch(std::vector<std::vector<std::string>>& sequenceNames,
                   std::vector<std::vector<std::string>>& sequences,
                   std::vector<std::vector<float>>& abundances,
                   Rcpp::List options);

    Rcpp::List removeChimeras();
    Rcpp::List removeChimeras(std::vector<std::string>& refNames,
                              std::vector<std::string>& refSequences);

private:

    Vsearch_Options* opts;
    bool dereplicate, bySample;
    int processors;

    // for reference and denovo without samples the size of pNames, pSequences,
    // and pAbundances will be 1
    std::vector<std::vector<std::string>>* pNames;
    std::vector<std::vector<std::string>>* pSequences;
    std::vector<std::vector<float>>* pAbundances;
    std::vector<std::string>* pRefNames;
    std::vector<std::string>* pRefSequences;

    void sortDescending(std::vector<std::string>& sequenceNames,
                        std::vector<std::string>& sequences,
                        std::vector<float>& abunds);

    bool contains(std::string s, Rcpp::CharacterVector& nv);

    Rcpp::List createReferencesProcesses();
    Rcpp::List createVsearchResults(std::vector<ChimeHit2>);
};
/******************************************************************************/
#endif // CHIMERAVSEARCH_H
