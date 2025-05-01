#include <testthat.h>
#include "kmer.h"

context("Kmer class C++ unit tests") {

    test_that("test getKmerNumber") {

        Kmer kmer(8);

        string seq = "TACGGAGGAGTGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGNTCGTTAA";

        int numKmers = seq.length() - 8 + 1;

        expect_true(kmer.getKmerNumber(seq, 0) == 50826);
        expect_true(kmer.getKmerNumber(seq, 1) == 6696);
        expect_true(kmer.getKmerNumber("NTCGTTAA", 0) == 14064);
    }
}
