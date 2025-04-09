#include <testthat.h>
#include "path.h"

context("PathData class C++ unit tests") {

    test_that("tests") {

        PathData path;

        expect_true(path.Start == 0);
        expect_true(path.Front == 0);
        expect_true(path.Back == 0);
        expect_true(path.Bytes == 0);

        path.Alloc(1000);

        expect_true(path.Bytes == 1001);
        expect_true(path.Start == 0);
        expect_false(path.Front == 0);
        expect_false(path.Back == 0);

        path.Alloc(1500);

        expect_true(path.Bytes == 1501);
        expect_true(path.Start == 0);
        expect_false(path.Front == 0);
        expect_false(path.Back == 0);

        path.Free();

        expect_true(path.Start == 0);
        expect_true(path.Front == 0);
        expect_true(path.Back == 0);
    }
}
