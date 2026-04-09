# Get path to rchime example files

rchime comes bundled with some example files in its `inst/extdata`
directory. This function make them easy to access.

## Usage

``` r
rchime_example(file = NULL)
```

## Arguments

- file:

  Name of file. If `NULL`, the example files will be listed.

## Examples

``` r
rchime_example()
#> [1] "/home/runner/work/_temp/Library/rchime/extdata"
rchime_example("reference.rds")
#> [1] "/home/runner/work/_temp/Library/rchime/extdata/reference.rds"
```
