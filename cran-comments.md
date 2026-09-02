## Test environments

- local OS package environment: macOS Tahoe, R 4.6.0
- GitHub Actions:
  - R CMD check of Ubuntu-latest (R-release, R-devel), Windows-latest, macOS-latest, macOS-14.
  - Test Coverage
  - Linting
  - Pkgdown build and deployment
  - Dependency only build and check

## R CMD check results

0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔

## Notes for Reviewer

All automated checks on GitHub Actions pass cleanly.

## Submission - 0.1.3

- Adds rchime::strollur_multi_sample_small() function to create the [strollur object](https://mothur.org/strollur/reference/strollur.html) used in examples.
- Adds internal functions to build the [strollur objects](https://mothur.org/strollur/reference/strollur.html) used in testing.
- Removes prebuilt [strollur objects](https://mothur.org/strollur/reference/strollur.html) from inst/extdata.


## Submission - 0.1.2

* Updates version from 0.1.1 to 0.1.2.
* Removes configure, cleanup and makevars files.
* Corrects compilation error for r-release-macos-x86_64, r-oldrel-macos-arm64, r-oldrel-macos-x86_64
* Adds old mac releases and extra flags to R-CMD-check github actions for better CRAN error detection.
* Updates example files for strollur v0.1.3 compatibility. 

## Submission - 0.1.1

* Updates version from 0.1.0 to 0.1.1.
* Corrects address sanitizing issue found on m1-san.
* Updates `rchime()` function to return a [strollur object](https://mothur.org/strollur/reference/strollur.html).

## Resubmission - 0.1.0

- Updates description field adding details about implementation and including 
    references to the work of Rognes et al. (2016) and Edgar,R.C. et al. (2011)
- Adds return value to rchime_example function.   

## Resubmission - 0.1.0

- Updates date field.
    - The Date field is over a month old.
- Corrects misspelled denovo to de novo.
    - Possibly misspelled words in DESCRIPTION:
        - Amplicon (3:40, 17:51)
        - denovo (17:103)

## Initial Submission - 0.1.0

Intial CRAN submission.
