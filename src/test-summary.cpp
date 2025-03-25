// Unit tests for Summary class

#include <testthat.h>
#include "summary.h"

context("Summary class C++ unit tests") {

    test_that("test summarizeFasta") {

        vector<int> dlengths = {250, 275, 233, 24, 240, 275};
        vector<int> dstarts = {1, 2, 1, 5, 50, 3};
        vector<int> dends = {213, 285, 243, 50, 300, 279};
        vector<int> dambigs = {1, 0, 0, 5, 2, 1};
        vector<int> dpolymers = {2, 1, 2, 2, 3, 7};
        vector<int> dnumns = {0, 0, 1, 0, 0, 2};

        Rcpp::DataFrame summary_report = DataFrame::create(
            Named("lengths") = dlengths,
            _["starts"] = dstarts,
            _["ends"] = dends,
            _["ambigs"] = dambigs,
            _["polymers"] = dpolymers,
            _["numns"] = dnumns);

        Rcpp::IntegerVector counts;
        counts.push_back(1);
        counts.push_back(10);
        counts.push_back(15);
        counts.push_back(1);
        counts.push_back(20);
        counts.push_back(1);

        // single thread
        Summary* summary = new Summary(1);

        Rcpp::DataFrame results = summary->summarizeFasta(summary_report,
                                                          counts);
        Rcpp::NumericVector starts = results["starts"];
        Rcpp::NumericVector ends = results["ends"];
        Rcpp::NumericVector ambigs = results["ambigs"];
        Rcpp::NumericVector nbases = results["nbases"];
        Rcpp::NumericVector polymers = results["polymers"];
        Rcpp::StringVector numseqs = results["numseqs"];

        // minimum
        expect_true(starts[0] == 1);
        // maximum
        expect_true(starts[6] == 50);

        // minimum
        expect_true(ends[0] == 50);
        // maximum
        expect_true(ends[6] == 300);

        // minimum
        expect_true(ambigs[0] == 0);
        // maximum
        expect_true(ambigs[6] == 5);

        // 25%
        expect_true(nbases[2] == 233);
        // 75%
        expect_true(nbases[4] == 250);

        // 25%
        expect_true(polymers[2] == 2);
        // 75%
        expect_true(polymers[4] == 3);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;

        // multi thread
        summary = new Summary(10);

        results = summary->summarizeFasta(summary_report, counts);
        starts = results["starts"];
        ends = results["ends"];
        ambigs = results["ambigs"];
        nbases = results["nbases"];
        polymers = results["polymers"];
        numseqs = results["numseqs"];

        // minimum
        expect_true(starts[0] == 1);
        // maximum
        expect_true(starts[6] == 50);

        // minimum
        expect_true(ends[0] == 50);
        // maximum
        expect_true(ends[6] == 300);

        // 25%
        expect_true(nbases[2] == 233);
        // 75%
        expect_true(nbases[4] == 250);

        // 25%
        expect_true(polymers[2] == 2);
        // 75%
        expect_true(polymers[4] == 3);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;
    }
}
