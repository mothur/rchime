#include <Rcpp.h>
#include "summary.h"

//' @name summarize_fasta
//' @title summarize_fasta
//' @rdname summarize_fasta
//' @param summary_report DataFrame, containing fasta summary data
//' @param count NumericalVector, containing abundances of sequences
//' @param processors Integer, number of cores to use. Default = all available
//' @description Summarizes a fasta dataset
//' @return DataFrame with summary values
//[[Rcpp::export]]
Rcpp::DataFrame summarize_fasta(Rcpp::DataFrame& summary_report,
                               Rcpp::IntegerVector& count, int processors) {

    // create Summary object
    Summary* summary = new Summary(processors);

    Rcpp::DataFrame results = summary->summarizeFasta(summary_report, count);

    delete summary;

    return results;
}
