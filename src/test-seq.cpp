#include <testthat.h>
#include "seq.h"

context("SeqData class C++ unit tests") {

    test_that("test constructors") {

        SeqData seq;
        expect_true(seq.getName() == "");
        expect_true(seq.getSeq() == "");
        expect_true(seq.getIndex() == UINT_MAX);
        expect_true(seq.getAbund() == 0);
        
        SeqData seq2("seq1", "ATGCCGGTTA", 15, 0);

        expect_true(seq2.getName() == "seq1");
        expect_true(seq2.getSeq() == "ATGCCGGTTA");
        expect_true(seq2.getIndex() == 0);
        expect_true(seq2.getAbund() == 15);

        seq2.setSeq("GGCGCCTTTAN");
        expect_true(seq2.getSeq() == "GGCGCCTTTAN");
        expect_true(seq2.getSeqLength() == 11);
    }
}
