// Rcpp functions using vector's as inputs

#include "chimera_vsearch.h"

// ============================================================================
//' @title Detects chimeras using a reference based approach
//' @name rchimeReference
//' @rdname rchimeReference
//'
//' @description detects chimeras in your data using a reference based approach.
//'
//' @param sequence_name, vector of strings containing sequence names
//' @param sequence, vector of strings containing sequence nucleotide data
//' @param reference_name, vector of strings containing reference sequences names
//' @param reference_sequence, vector of strings containing reference sequences
//' @param options list containing parameter options. Default = NULL.
//'
//' @examples
//'
//' # Read in data.frames containing sequence data
//'
//' fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
//' reference_data <- readRDS(rchime_example("reference.rds"))
//'
//' # Detect chimeras in first 100 sequences
//'
//' results <- rchimeReference(sequence_name = fasta_data$sequence_name[1:100],
//'                            sequence = fasta_data$sequence[1:100],
//'                            reference_name = reference_data$sequence_name,
//'                            reference_sequence = reference_data$sequence)
//'
//' @seealso [rchime()]
//' @export
//[[Rcpp::export]]
Rcpp::List rchimeReference(std::vector<std::string> sequence_name,
                           std::vector<std::string> sequence,
                           std::vector<std::string> reference_name,
                           std::vector<std::string> reference_sequence,
                           Rcpp::Nullable<Rcpp::List> options = R_NilValue) {

    // check vector lengths match
    if (sequence_name.size() != sequence.size()) {
        std::string message = "sequence_name must be the same length as sequence.";
        message += " You provided " + toString(sequence_name.size());
        message += " sequence names, and " + toString(sequence.size());
        message += " sequences.";
        throw Rcpp::exception(message.c_str());
    }

    if (reference_name.size() != reference_sequence.size()) {
        std::string message = "reference_name must be the same length as reference_sequence.";
        message += " You provided " + toString(reference_name.size());
        message += " reference names, and " + toString(reference_sequence.size());
        message += " reference sequences.";
        throw Rcpp::exception(message.c_str());
    }

    // Use Rcpp::Environment and Rcpp::Function to call R code from C++.
    Rcpp::Environment parallelly_env = Rcpp::Environment::namespace_env("parallelly");
    Rcpp::Function availableCores = parallelly_env["availableCores"];

    Rcpp::List optionsReference = Rcpp::List::create();
    if (options.isNotNull()) {

        Rcpp::List opts_list = Rcpp::as<Rcpp::List>(options);

        // add required processors options if needed
        if (!opts_list.containsElementNamed("processors")) {
            std::vector<std::string> newNames = opts_list.attr("names");
            newNames.push_back("processors");
            opts_list.push_back(Rcpp::as<int>(availableCores()));
            opts_list.attr("names") = newNames;
        }

        // add required dereplicate options if needed
        if (!opts_list.containsElementNamed("dereplicate")) {
            std::vector<std::string> newNames = opts_list.attr("names");
            newNames.push_back("dereplicate");
            opts_list.push_back(false);
            opts_list.attr("names") = newNames;
        }
        optionsReference = opts_list;
    }else{
        std::vector<std::string> optionNames = {"processors", "dereplicate"};
        optionsReference.push_back(Rcpp::as<int>(availableCores()));
        optionsReference.push_back(false);
        optionsReference.attr("names") = optionNames;
    }

    std::vector<std::vector<std::string>> seqNames;
    seqNames.push_back(sequence_name);
    std::vector<std::vector<std::string>> seqs;
    seqs.push_back(sequence);
    std::vector<std::vector<float>> abunds;
    std::vector<float> abundances(sequence_name.size(), 1);
    abunds.push_back(abundances);

    ChimeraVsearch* chimera = new ChimeraVsearch(seqNames, seqs, abunds,
                                                 optionsReference);

    // results contains chimera_report (data.frame) and chimeras (vector of names)
    Rcpp::List results = chimera->removeChimeras(reference_name,
                                                 reference_sequence);

    delete chimera;

    return results;
}
// ============================================================================
//' @title Detects chimeras using a denovo approach
//' @name rchimeDenovoSingleSample
//' @rdname rchimeDenovoSingleSample
//'
//' @description detects chimeras in your data using a denovo approach.
//'
//' @param sequence_name, vector of strings containing sequence names
//' @param sequence, vector of strings containing sequence nucleotide data
//' @param abundance, vector of sequence abundances
//' @param options list containing parameter options. Default = NULL.
//'
//' @examples
//'
//' # Read in data.frames containing sequence and abundance data
//'
//' fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
//' abundance_data <- readRDS(rchime_example("single_sample_abundance.rds"))
//'
//' # Detect chimeras in first 100 sequences
//'
//' results <- rchimeDenovoSingleSample(sequence_name = fasta_data$sequence_name[1:100],
//'                                     sequence = fasta_data$sequence[1:100],
//'                                     abundance = abundance_data$abundance[1:100])
//'
//' @seealso [rchime()]
//' @export
//[[Rcpp::export]]
Rcpp::List rchimeDenovoSingleSample(std::vector<std::string> sequence_name,
                        std::vector<std::string> sequence,
                        std::vector<float> abundance,
                        Rcpp::Nullable<Rcpp::List> options = R_NilValue) {

    // check vector lengths match
    if (sequence_name.size() != sequence.size()) {
        std::string message = "sequence_name must be the same length as sequence.";
        message += " You provided " + toString(sequence_name.size());
        message += " sequence names, and " + toString(sequence.size());
        message += " sequences.";
        throw Rcpp::exception(message.c_str());
    }

    // check vector lengths match
    if (sequence_name.size() != abundance.size()) {
        std::string message = "sequence_name must be the same length as abundance.";
        message += " You provided " + toString(sequence_name.size());
        message += " sequence names, and " + toString(abundance.size());
        message += " abundances.";
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::List optionsDenovo = Rcpp::List::create();
    if (options.isNotNull()) {

        Rcpp::List opts_list = Rcpp::as<Rcpp::List>(options);

        // add required processors options if needed
        if (!opts_list.containsElementNamed("processors")) {
            std::vector<std::string> newNames = opts_list.attr("names");
            newNames.push_back("processors");
            opts_list.push_back(1);
            opts_list.attr("names") = newNames;
        }

        // add required dereplicate options if needed
        if (!opts_list.containsElementNamed("dereplicate")) {
            std::vector<std::string> newNames = opts_list.attr("names");
            newNames.push_back("dereplicate");
            opts_list.push_back(true);
            opts_list.attr("names") = newNames;
        }
        optionsDenovo = opts_list;
    }else{
        std::vector<std::string> optionNames = {"processors", "dereplicate"};
        optionsDenovo.push_back(1);
        optionsDenovo.push_back(true);
        optionsDenovo.attr("names") = optionNames;
    }

    std::vector<std::vector<std::string>> seqNames;
    seqNames.push_back(sequence_name);
    std::vector<std::vector<std::string>> seqs;
    seqs.push_back(sequence);
    std::vector<std::vector<float>> abunds;
    abunds.push_back(abundance);

    ChimeraVsearch* chimera = new ChimeraVsearch(seqNames, seqs, abunds,
                                                 optionsDenovo);

    Rcpp::List results = chimera->removeChimeras();

    delete chimera;

    return results;
}
// ============================================================================
//' @title Detects chimeras using a denovo method processing by sample.
//' @name rchimeDenovo
//' @rdname rchimeDenovo
//'
//' @param sequence_name, 2D vector of strings containing sequence names parsed by sample
//' @param sequence, 2D vector of strings containing sequence nucleotide data  parsed by sample
//' @param abundance, 2D vector of sequence abundances parsed by sample
//' @param options list containing parameter options. Default = NULL.
//'
//' @examples
//'
//' # Read in data.frames containing sequence and abundance data
//'
//' sequence_name <- readRDS(rchime_example("miseq_names_by_sample.rds"))
//' sequence <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
//' abundance <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))
//'
//' # Detect chimeras in first sample
//'
//' results <- rchimeDenovo(sequence_name = sequence_name[1],
//'                         sequence = sequence[1],
//'                         abundance = abundance[1])
//'
//' @seealso [rchime()]
//' @description detects chimeras from your data using a denovo method processing by sample.
//' @returns list()
//' @export
//[[Rcpp::export]]
Rcpp::List rchimeDenovo(std::vector<std::vector<std::string>> sequence_name,
                         std::vector<std::vector<std::string>> sequence,
                         std::vector<std::vector<float>> abundance,
                         Rcpp::Nullable<Rcpp::List> options = R_NilValue) {

    // check vector lengths match
    if (sequence_name.size() != sequence.size()) {
        std::string message = "sequence_name must be the same length as sequence.";
        throw Rcpp::exception(message.c_str());
    }

    // check vector lengths match
    if (sequence_name.size() != abundance.size()) {
        std::string message = "sequence_name must be the same length as abundance.";
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::List optionsDenovo = Rcpp::List::create();
    if (options.isNotNull()) {

        Rcpp::List opts_list = Rcpp::as<Rcpp::List>(options);

        // add required processors options if needed
        if (!opts_list.containsElementNamed("processors")) {
            std::vector<std::string> newNames = opts_list.attr("names");
            newNames.push_back("processors");
            opts_list.push_back(1);
            opts_list.attr("names") = newNames;
        }

        // add required dereplicate options if needed
        if (!opts_list.containsElementNamed("dereplicate")) {
            std::vector<std::string> newNames = opts_list.attr("names");
            newNames.push_back("dereplicate");
            opts_list.push_back(true);
            opts_list.attr("names") = newNames;
        }
        optionsDenovo = opts_list;
    }else{
        std::vector<std::string> optionNames = {"processors", "dereplicate"};
        optionsDenovo.push_back(1);
        optionsDenovo.push_back(true);
        optionsDenovo.attr("names") = optionNames;
    }

     ChimeraVsearch* chimera = new ChimeraVsearch(sequence_name,
                                                  sequence,
                                                  abundance,
                                                  optionsDenovo);

     Rcpp::List results = chimera->removeChimeras();

     delete chimera;

     return results;
}
// ============================================================================
