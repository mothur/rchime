// Unit tests for SequenceParser class

#include <testthat.h>
#include "sequenceparser.h"

context("SequenceParser class C++ unit tests") {

    test_that("test constructor") {

        Rcpp::Environment rchime_env("package:rchime");
        Environment dataset_env = rchime_env["sequence_data_vector"];
        Function newFastaInstance = dataset_env["new"];
        Rcpp::Environment fastaDatabase = newFastaInstance();

        vector<string> names, seqs, groups;

        names.push_back("seq1");
        seqs.push_back("ATTGCCGTA");
        groups.push_back("sample1");

        names.push_back("seq2");
        seqs.push_back("ATTGCCGCA");
        groups.push_back("sample2");

        names.push_back("seq3");
        seqs.push_back("ATATCCGTA");
        groups.push_back("sample1");

        names.push_back("seq4");
        seqs.push_back("ATTGTTGTA");
        groups.push_back("sample2");

        Function addSeqs = fastaDatabase["add_seqs"];
        addSeqs(names, seqs);

        Function addGroups = fastaDatabase["set_group_assignments"];
        addGroups(names, groups);

        SequenceParser parser(fastaDatabase, true);

        vector<vector<string>> out_names = parser.getNames(nullVector);

        expect_true(names.size() == out_names[0].size());
        expect_true(names[0] == out_names[0][0]);
        expect_true(names[1] == out_names[0][1]);
        expect_true(names[2] == out_names[0][2]);
        expect_true(names[3] == out_names[0][3]);
    }
}
