// Rcpp functions using vector's as inputs

#include "chimera_vsearch.h"

// ============================================================================
//' @title Detects chimeras using a reference based approach
//' @name rchimeReference
//' @rdname rchimeReference
//'
//' @description detects chimeras in your data using a reference based approach.
//'
//' @param sequence_names, vector of strings containing sequence names
//' @param sequences, vector of strings containing sequence nucleotide data
//' @param reference_names, vector of strings containing reference sequences names
//' @param reference_sequences, vector of strings containing reference sequences
//' @param options list containing parameter options. Default = NULL.
//'
//' @examples
//'
//' # Read in data.frames containing sequence data
//'
//' fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
//' reference_data <- readRDS(rchime_example("reference.rds"))
//'
//' # Detect chimeras
//'
//' results <- rchimeReference(sequence_names = fasta_data$sequence_names,
//'                            sequences = fasta_data$sequences,
//'                            reference_names = reference_data$sequence_names,
//'                            reference_sequences = reference_data$sequences)
//'
//' @seealso [rchime()]
//' @export
//[[Rcpp::export]]
Rcpp::List rchimeReference(std::vector<std::string> sequence_names,
                           std::vector<std::string> sequences,
                           std::vector<std::string> reference_names,
                           std::vector<std::string> reference_sequences,
                           Rcpp::Nullable<Rcpp::List> options = R_NilValue) {

    // check vector lengths match
    if (sequence_names.size() != sequences.size()) {
        std::string message = "sequence_names must be the same length as sequences.";
        message += " You provided " + toString(sequence_names.size());
        message += " sequence names, and " + toString(sequences.size());
        message += " sequences.";
        throw Rcpp::exception(message.c_str());
    }

    if (reference_names.size() != reference_sequences.size()) {
        std::string message = "reference_names must be the same length as reference_sequences.";
        message += " You provided " + toString(reference_names.size());
        message += " reference names, and " + toString(reference_sequences.size());
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
    seqNames.push_back(sequence_names);
    std::vector<std::vector<std::string>> seqs;
    seqs.push_back(sequences);
    std::vector<std::vector<float>> abunds;
    std::vector<float> abundances(sequence_names.size(), 1);
    abunds.push_back(abundances);

    ChimeraVsearch* chimera = new ChimeraVsearch(seqNames, seqs, abunds,
                                                 optionsReference);

    // results contains chimera_report (data.frame) and chimeras (vector of names)
    Rcpp::List results = chimera->removeChimeras(reference_names,
                                                 reference_sequences);

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
//' @param sequence_names, vector of strings containing sequence names
//' @param sequences, vector of strings containing sequence nucleotide data
//' @param abundances, vector of sequence abundances
//' @param options list containing parameter options. Default = NULL.
//'
//' @examples
//'
//' # Read in data.frames containing sequence and abundance data
//'
//' fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
//' abundance_data <- readRDS(rchime_example("single_sample_abundance.rds"))
//'
//' # Detect chimeras
//'
//' results <- rchimeDenovoSingleSample(sequence_names = fasta_data$sequence_names,
//'                                     sequences = fasta_data$sequences,
//'                                     abundances = abundance_data$abundances)
//'
//' @seealso [rchime()]
//' @export
//[[Rcpp::export]]
Rcpp::List rchimeDenovoSingleSample(std::vector<std::string> sequence_names,
                        std::vector<std::string> sequences,
                        std::vector<float> abundances,
                        Rcpp::Nullable<Rcpp::List> options = R_NilValue) {

    // check vector lengths match
    if (sequence_names.size() != sequences.size()) {
        std::string message = "sequence_names must be the same length as sequences.";
        message += " You provided " + toString(sequence_names.size());
        message += " sequence names, and " + toString(sequences.size());
        message += " sequences.";
        throw Rcpp::exception(message.c_str());
    }

    // check vector lengths match
    if (sequence_names.size() != abundances.size()) {
        std::string message = "sequence_names must be the same length as abundances.";
        message += " You provided " + toString(sequence_names.size());
        message += " sequence names, and " + toString(abundances.size());
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
            opts_list.push_back(false);
            opts_list.attr("names") = newNames;
        }
        optionsDenovo = opts_list;
    }else{
        std::vector<std::string> optionNames = {"processors", "dereplicate"};
        optionsDenovo.push_back(1);
        optionsDenovo.push_back(false);
        optionsDenovo.attr("names") = optionNames;
    }

    std::vector<std::vector<std::string>> seqNames;
    seqNames.push_back(sequence_names);
    std::vector<std::vector<std::string>> seqs;
    seqs.push_back(sequences);
    std::vector<std::vector<float>> abunds;
    abunds.push_back(abundances);

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
//' @param sequence_names, 2D vector of strings containing sequence names parsed by sample
//' @param sequences, 2D vector of strings containing sequence nucleotide data  parsed by sample
//' @param abundances, 2D vector of sequence abundances parsed by sample
//' @param options list containing parameter options. Default = NULL.
//'
//' @examples
//'
//' # Read in data.frames containing sequence and abundance data
//'
//' sequence_names <- readRDS(rchime_example("miseq_names_by_sample.rds"))
//' sequences <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
//' abundances <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))
//'
//' # Detect chimeras
//'
//' options <- rchime_options(dereplicate = TRUE)
//' results <- rchimeDenovo(sequence_names = sequence_names,
//'                         sequences = sequences,
//'                         abundances = abundances,
//'                         options)
//'
//' @seealso [rchime()]
//' @description detects chimeras from your data using a denovo method processing by sample.
//' @returns list()
//' @export
//[[Rcpp::export]]
Rcpp::List rchimeDenovo(std::vector<std::vector<std::string>> sequence_names,
                         std::vector<std::vector<std::string>> sequences,
                         std::vector<std::vector<float>> abundances,
                         Rcpp::Nullable<Rcpp::List> options = R_NilValue) {

    // check vector lengths match
    if (sequence_names.size() != sequences.size()) {
        std::string message = "sequence_names must be the same length as sequences.";
        throw Rcpp::exception(message.c_str());
    }

    // check vector lengths match
    if (sequence_names.size() != abundances.size()) {
        std::string message = "sequence_names must be the same length as abundances.";
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
            opts_list.push_back(false);
            opts_list.attr("names") = newNames;
        }
        optionsDenovo = opts_list;
    }else{
        std::vector<std::string> optionNames = {"processors", "dereplicate"};
        optionsDenovo.push_back(1);
        optionsDenovo.push_back(false);
        optionsDenovo.attr("names") = optionNames;
    }

     ChimeraVsearch* chimera = new ChimeraVsearch(sequence_names,
                                                  sequences,
                                                  abundances,
                                                  optionsDenovo);

     Rcpp::List results = chimera->removeChimeras();

     delete chimera;

     return results;
}
// ============================================================================
