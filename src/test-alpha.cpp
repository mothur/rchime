#include <testthat.h>
#include "alpha.h"

context("Alpha C++ unit tests") {

    test_that("test g_CharToLetterNucleo") {

        for (int i = 0; i < 65; i++) {
            expect_true(g_CharToLetterNucleo[i] == INVALID_LETTER);
        }
        for (int i = 118; i < 256; i++) {
            expect_true(g_CharToLetterNucleo[i] == INVALID_LETTER);
        }

        // A / a
        expect_true(g_CharToLetterNucleo[65] == 0);
        expect_true(g_CharToLetterNucleo[97] == 0);

        // C / c
        expect_true(g_CharToLetterNucleo[67] == 1);
        expect_true(g_CharToLetterNucleo[99] == 1);

        // G / g
        expect_true(g_CharToLetterNucleo[71] == 2);
        expect_true(g_CharToLetterNucleo[103] == 2);

        // T & U / t & u
        expect_true(g_CharToLetterNucleo[84] == 3);
        expect_true(g_CharToLetterNucleo[85] == 3);
        expect_true(g_CharToLetterNucleo[116] == 3);
        expect_true(g_CharToLetterNucleo[117] == 3);
    }

    test_that("test g_IsACGTU") {

        for (int i = 0; i < 65; i++) {
            expect_false(g_IsACGTU[i]);
        }
        for (int i = 118; i < 256; i++) {
            expect_false(g_IsACGTU[i]);
        }

        // A / a
        expect_true(g_IsACGTU[65]);
        expect_true(g_IsACGTU[97]);

        // C / c
        expect_true(g_IsACGTU[67]);
        expect_true(g_IsACGTU[99]);

        // G / g
        expect_true(g_IsACGTU[71]);
        expect_true(g_IsACGTU[103]);

        // T & U / t & u
        expect_true(g_IsACGTU[84]);
        expect_true(g_IsACGTU[85]);
        expect_true(g_IsACGTU[116]);
        expect_true(g_IsACGTU[117]);
    }
}
