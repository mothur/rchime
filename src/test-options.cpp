#include <testthat.h>
#include "options.h"

context("Options class C++ unit tests") {

    test_that("test defaults") {

        Options opts;

        expect_true(opts.getAbskew() == 2.0);
        expect_false(opts.getChimealns());
        expect_true(opts.getChunks() == 4);
        expect_true(opts.getDn() == 1.4);
        expect_true(opts.getIdsmoothwindow() == 32);
        expect_true(opts.getMinchunk() == 64);

        expect_true(opts.getMindiffs() == 3);
        expect_true(opts.getSkipgaps());
        expect_true(opts.getSkipgaps2());
        expect_true(opts.getMindiv() == 0.5);
        expect_true(opts.getMinh() == 0.3);
        expect_true(opts.getMinlen() == 10);
        expect_true(opts.getMaxlen() == 10000);

        expect_true(opts.getMaxp() == 2);
        expect_true(opts.getQueryfract() == 0.5);
        expect_true(opts.getXa() == 1.0);
        expect_true(opts.getXn() == 8.0);
    }
}
