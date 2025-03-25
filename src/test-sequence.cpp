// Unit tests for Sequence class

#include <testthat.h>
#include "sequence.h"

context("Sequence class C++ unit tests") {

    test_that("test constructors") {
        Sequence seq;
        expect_true(seq.getName() == "");
        expect_true(seq.getUnaligned() == "");
        expect_true(seq.getAligned() == "");
        expect_true(seq.getPairwise() == "");
        expect_true(seq.getComment() == "");
        expect_true(seq.getNumBases() == 0);
        expect_true(seq.getAlignLength() == 0);
        expect_true(seq.getStartPos() == 1);
        expect_true(seq.getEndPos() == 0);
        expect_true(seq.getLongHomoPolymer() == 1);
        expect_true(seq.getAmbigBases() == 0);

        Sequence seq2("testSeq", "A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-T-");
        expect_true(seq2.getName() == "testSeq");
        expect_true(seq2.getUnaligned() == "AGCGAGAAGTATGCGGAATGCGT");
        expect_true(seq2.getAligned() ==
            "A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-T-");
        expect_true(seq2.getPairwise() == "");
        expect_true(seq2.getComment() == "");
        expect_true(seq2.getNumBases() == 23);
        expect_true(seq2.getAlignLength() == 43);
        expect_true(seq2.getStartPos() == 1);
        expect_true(seq2.getEndPos() == 42);
        expect_true(seq2.getLongHomoPolymer() == 2);
        expect_true(seq2.getAmbigBases() == 0);
    }

    test_that("test getters and setters") {
        Sequence seq;

        seq.setName("testSeq");
        expect_true(seq.getName() == "testSeq");

        seq.setAligned("A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-N-");
        expect_true(seq.getAligned() ==
            "A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-N-");
        expect_true(seq.getUnaligned() == "AGCGAGAAGTATGCGGAATGCGN");
        seq.setUnaligned("AGCGAGAAGTATGCGGAATGCGN");
        expect_true(seq.getUnaligned() == "AGCGAGAAGTATGCGGAATGCGN");

        seq.setPairwise("ATGCGC");
        expect_true(seq.getPairwise() == "ATGCGC");

        seq.setComment("this is my really cool sequence");
        expect_true(seq.getComment() == "this is my really cool sequence");

        expect_true(seq.getNumBases() == 23);
        expect_true(seq.getAlignLength() == 43);
        expect_true(seq.getStartPos() == 1);
        expect_true(seq.getEndPos() == 42);
        expect_true(seq.getLongHomoPolymer() == 2);
        expect_true(seq.getAmbigBases() == 1);
        expect_true(seq.getNumNs() == 1);
        expect_true(seq.isAligned());

        expect_true(seq.getInlineSeq() ==
            "testSeq\tA-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-N-");

        seq.setAligned("-A-GC--GA-G-A-A-G-U-A--TG-C--GG-A-ATG-C-G-N-");
        expect_true(seq.getAligned() ==
            ".A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-N.");
        expect_true(seq.getStartPos() == 2);
    }

    test_that("test pads and filters") {
        Sequence seq("testSeq", "A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-T-");

        seq.padToPos(10);
        expect_true(seq.getAligned() ==
            ".........G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-T-");

        seq.padFromPos(40);
        expect_true(seq.getAligned() ==
            ".........G-A-A-G-T-A--TG-C--GG-A-ATG-C-G..-");


        seq.filterToPos(15);
        expect_true(seq.getAligned() ==
            "...............G-T-A--TG-C--GG-A-ATG-C-G..-");

        seq.filterFromPos(35);
        expect_true(seq.getAligned() ==
            "...............G-T-A--TG-C--GG-A-AT........");

        Sequence seq2("testSeq", "....---A-G-T-A--TG-C--GG-A-ATG-C-G-T---...");
        seq2.filterToPos(15);
        expect_true(seq2.getAligned() ==
            "................TG-C--GG-A-ATG-C-G-T---...");

        seq2.filterFromPos(25);
        expect_true(seq2.getAligned() ==
            "................TG-C--GG..................");
    }

    test_that("test class functions") {
        Sequence seq("testSeq", "A-GC--GA-G-A-A-G-T-A--TG-C--GG-A-ATG-C-G-N-");

        seq.reverseComplement();
        expect_true(seq.getAligned() ==
            "NCGCATTCCGCATACTTCTCGCT");

        seq.trim(10);
        expect_true(seq.getAligned() ==
            "NCGCATTCCG");

        seq.removeAmbigBases();
        expect_true(seq.getAligned() ==
            "-CGCATTCCG");

        expect_true(seq.convert2ints() ==
            "121033112");

        Sequence seq2("testSeq", "NCGCAUUCWG");
        expect_true(seq2.convert2ints() == "4121033142");
        expect_true(seq2.isAligned() == false);
    }
}
