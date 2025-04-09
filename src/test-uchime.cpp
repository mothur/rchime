#include <testthat.h>
#include "uchime.h"

context("Uchime global functions C++ unit tests") {

    test_that("test pom tof yon") {

        expect_true(yon(true) == 'Y');
        expect_true(yon(false) == 'N');

        expect_true(tof(true) == 'T');
        expect_true(tof(false) == 'F');

        expect_true(pom(true) == '+');
        expect_true(pom(false) == '-');
    }

    test_that("test isacgt isgap") {

        expect_true(isacgt('A'));
        expect_true(isacgt('T'));
        expect_true(isacgt('G'));
        expect_true(isacgt('C'));
        expect_false(isacgt('N'));

        expect_true(isgap('-'));
        expect_true(isgap('.'));
        expect_false(isgap('A'));
    }
}
