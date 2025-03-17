#include "chimera_uchime.h"
#include "chimera.h"

// ============================================================================
//' @title chimeraUchimeReference
//' @name chimeraUchimeReference
//' @rdname chimeraUchimeReference
//' @param dataset R6 'sequence_dataset' object containing sequence data
//' @param reference R6 'sequence_dataset' object containing sequence data
//' @param options list containing parameter options
//' @seealso [chimera_uchime()]
//' @description detects and removes chimeras from your data using a reference
//'  dataset.
//[[Rcpp::export]]
Rcpp::List chimeraUchimeReference(Rcpp::Environment& dataset,
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
//' @title chimeraUchime
//' @name chimeraUchime
//' @rdname chimeraUchime
//' @param dataset R6 'sequence_dataset' object containing sequence data
//' @param options list containing parameter options
//' @seealso [chimera_uchime()]
//' @description detects and removes chimeras from your data using a denovo
//' approach.
//[[Rcpp::export]]
Rcpp::List chimeraUchime(Rcpp::Environment& dataset,
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
