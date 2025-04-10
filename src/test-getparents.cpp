#include <testthat.h>
#include "getparents.h"

context("GetParents class C++ unit tests") {

    test_that("test getChunkInfo") {

        Options opt;
        GetParents parentFinder(&opt);

        unsigned ChunkLength;
	    vector<unsigned> ChunkLos = parentFinder.getChunkInfo(252, ChunkLength);

        expect_true(ChunkLength == 64);
        expect_true(ChunkLos[0] == 0);
        expect_true(ChunkLos[1] == 64);
        expect_true(ChunkLos[2] == 128);
        expect_true(ChunkLos[3] == 187);

        ChunkLos = parentFinder.getChunkInfo(195, ChunkLength);

        expect_true(ChunkLength == 64);
        expect_true(ChunkLos[0] == 0);
        expect_true(ChunkLos[1] == 64);
        expect_true(ChunkLos[2] == 128);
        expect_true(ChunkLos[3] == 130);
    }

    test_that("test getWord") {
        Options opt;
        GetParents parentFinder(&opt);

        string sequence = "ATGGCATGTTACC";

        expect_true(parentFinder.getWord(sequence, 0) == 14926);
        expect_true(parentFinder.getWord(sequence, 1) == 59707);
        expect_true(parentFinder.getWord(sequence, 2) == 42223);
        expect_true(parentFinder.getWord(sequence, 3) == 37820);
        expect_true(parentFinder.getWord(sequence, 4) == 20209);
        expect_true(parentFinder.getWord(sequence, 5) == 15301);
        expect_true(parentFinder.getWord(sequence, 6) == 61204);
    }
}
