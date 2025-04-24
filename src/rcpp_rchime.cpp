#include "chimera_uchime.h"
#include "chimera.h"

// ============================================================================
//' @title rchimeReference
//' @name rchimeReference
//' @rdname rchimeReference
//' @param dataset R6 'sequence_dataset' object containing sequence data
//' @param reference R6 'sequence_dataset' object containing sequence data
//' @param options list containing parameter options
//' @seealso [rchime()]
//' @description detects and removes chimeras from your data using a reference
//'  dataset.
//[[Rcpp::export]]
Rcpp::List rchimeReference(Rcpp::Environment& dataset,
                       Rcpp::Environment& reference,
                       Rcpp::List options) {

     bool dereplicate = Rcpp::as<bool>(options["dereplicate"]);
     bool silent = Rcpp::as<bool>(options["silent"]);
     int processors = Rcpp::as<int>(options["processors"]);

     Rcpp::Function getNumGroups = dataset["get_num_groups"];

     int numGroups = Rcpp::as<int>(getNumGroups());
     bool hasGroupData = true;
     if (numGroups == 0) {
         hasGroupData = false;
     }

     Chimera* chimera = new ChimeraUchime(dereplicate, processors, silent,
                                            hasGroupData, options);

     Rcpp::List chimera_report = chimera->removeChimeras(dataset, reference);

     delete chimera;

     return chimera_report;
 }

// ============================================================================
//' @title rchimeDenovo
//' @name rchimeDenovo
//' @rdname rchimeDenovo
//' @param dataset R6 'sequence_dataset' object containing sequence data
//' @param options list containing parameter options
//' @seealso [rchime()]
//' @description detects and removes chimeras from your data using a denovo
//' approach.
//[[Rcpp::export]]
Rcpp::List rchimeDenovo(Rcpp::Environment& dataset,
                       Rcpp::List options) {

     bool dereplicate = Rcpp::as<bool>(options["dereplicate"]);
     bool silent = Rcpp::as<bool>(options["silent"]);
     int processors = Rcpp::as<int>(options["processors"]);

     Rcpp::Function getNumGroups = dataset["get_num_groups"];

     int numGroups = Rcpp::as<int>(getNumGroups());
     bool hasGroupData = true;
     if (numGroups == 0) {
         hasGroupData = false;
     }

     Chimera* chimera = new ChimeraUchime(dereplicate, processors, silent,
                                            hasGroupData, options);

     Rcpp::List chimera_report = chimera->removeChimeras(dataset);

     delete chimera;

     return chimera_report;
}
// ============================================================================
