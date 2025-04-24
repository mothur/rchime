#include <testthat.h>
#include "alnparams.h"

context("AlnParams class C++ unit tests") {

    test_that("test constructor") {

        AlnParams alignParameters;

        expect_true(alignParameters.OpenA == -10.0);
        expect_true(alignParameters.OpenB == -10.0);
        expect_true(alignParameters.ExtA == -1.0);
        expect_true(alignParameters.ExtB == -1.0);

        expect_true(alignParameters.LOpenA == -0.5);
        expect_true(alignParameters.LOpenB == -0.5);
        expect_true(alignParameters.ROpenA == -0.5);
        expect_true(alignParameters.ROpenB == -0.5);

        expect_true(alignParameters.LExtA == -0.5);
        expect_true(alignParameters.LExtB == -0.5);
        expect_true(alignParameters.RExtA == -0.5);
        expect_true(alignParameters.RExtB == -0.5);

        expect_true(alignParameters.SubstMx.matrix.size() == 272);
    }
}
