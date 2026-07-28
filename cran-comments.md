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

## Test environments

- local OS package environment: macOS Tahoe, R 4.6.0
- GitHub Actions:
  - R CMD check of Ubuntu-latest (R-release, R-devel), Windows-latest, macOS-latest.
  - Test Coverage
  - Linting
  - Pkgdown build and deployment
  - Dependency only build and check

## R CMD check results

0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔

## Notes for Reviewer

All automated checks on GitHub Actions pass cleanly.
