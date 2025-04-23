#include <testthat.h>
#include "mxmatrix.h"

context("MyMatrix class C++ unit tests") {

    test_that("test defaults") {

        MxFloatMatrix m(10, 20, 0);

        expect_true(m.matrix.size() == 36);

        m.resize(20, 30);

        expect_true(m.matrix.size() == 46);
    }
}
