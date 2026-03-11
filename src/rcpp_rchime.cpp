#include "chimera_vsearch.h"

// ============================================================================
//' @title rchimeReference
//' @name rchimeReference
//' @rdname rchimeReference
//' @param sequence_names, vector of strings containing sequence names
//' @param sequences, vector of strings containing sequence nucleotide data
//' @param abundances, vector of sequence abundances
//' @param reference_names, vector of strings containing reference sequences names
//' @param reference_sequences, vector of strings containing reference sequences
//' @param options list containing parameter options
//' @seealso [rchime()]
//' @description detects chimeras from your data using a reference dataset.
//[[Rcpp::export]]
Rcpp::List rchimeReference(std::vector<std::string> sequence_names,
                           std::vector<std::string> sequences,
                           std::vector<float> abundances,
                           std::vector<std::string> reference_names,
                           std::vector<std::string> reference_sequences,
                           Rcpp::List options) {

    std::vector<std::vector<std::string>> seqNames;
    seqNames.push_back(sequence_names);
    std::vector<std::vector<std::string>> seqs;
    seqs.push_back(sequences);
    std::vector<std::vector<float>> abunds;
    abunds.push_back(abundances);

    ChimeraVsearch* chimera = new ChimeraVsearch(seqNames, seqs, abunds,
                                                 options);

    // results contains chimera_report (data.frame) and chimeras (vector of names)
    Rcpp::List results = chimera->removeChimeras(reference_names,
                                                 reference_sequences);

    delete chimera;

    return results;
}
// ============================================================================
//' @title rchimeDenovoSingleSample
//' @name rchimeDenovoSingleSample
//' @rdname rchimeDenovoSingleSample
//' @param sequence_names, vector of strings containing sequence names
//' @param sequences, vector of strings containing sequence nucleotide data
//' @param abundances, vector of sequence abundances
//' @param options list containing parameter options
//' @seealso [rchime()]
//' @description detects chimeras from your data using a denovo method.
//[[Rcpp::export]]
Rcpp::List rchimeDenovoSingleSample(std::vector<std::string> sequence_names,
                        std::vector<std::string> sequences,
                        std::vector<float> abundances,
                        Rcpp::List options) {


    std::vector<std::vector<std::string>> seqNames;
    seqNames.push_back(sequence_names);
    std::vector<std::vector<std::string>> seqs;
    seqs.push_back(sequences);
    std::vector<std::vector<float>> abunds;
    abunds.push_back(abundances);

    ChimeraVsearch* chimera = new ChimeraVsearch(seqNames, seqs, abunds,
                                                 options);

    Rcpp::List results = chimera->removeChimeras();

    delete chimera;

    return results;
}
//' @title rchimeDenovo
//' @name rchimeDenovo
//' @rdname rchimeDenovo
//' @param sequence_names, 2D vector of strings containing sequence names parsed by sample
//' @param sequences, 2D vector of strings containing sequence nucleotide data  parsed by sample
//' @param abundances, 2D vector of sequence abundances parsed by sample
//' @param options list containing parameter options
//' @seealso [rchime()]
//' @description detects chimeras from your data using a denovo method processing by sample.
//[[Rcpp::export]]
Rcpp::List rchimeDenovo(std::vector<std::vector<std::string>> sequence_names,
                         std::vector<std::vector<std::string>> sequences,
                         std::vector<std::vector<float>> abundances,
                         Rcpp::List options) {


     ChimeraVsearch* chimera = new ChimeraVsearch(sequence_names,
                                                  sequences,
                                                  abundances,
                                                  options);

     Rcpp::List results = chimera->removeChimeras();

     delete chimera;

     return results;
}
// ============================================================================

