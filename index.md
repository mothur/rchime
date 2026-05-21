# rchime

[![pkgdown](https://github.com/SchlossLab/rchime/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/SchlossLab/rchime/actions/workflows/pkgdown.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/rchime)](https://CRAN.R-project.org/package=rchime)

## Overview

The *rchime* package allows you to detect and remove chimeras from your
dataset using a denovo approach or alternatively a reference model. This
package uses code from the
[vsearch](https://github.com/torognes/vsearch) tools.

- [`rchime()`](http://mothur.org/rchime/reference/rchime.md) detect and
  remove chimeras from your [strollur](https://mothur.org/strollur/)
  dataset object or data.frame

## Installation

You can install the CRAN version with:

``` r

install.packages("rchime")
```

## Development version

You can install the development version of rchime from
[GitHub](https://github.com/mothur/rchime) with:

``` r

# install.packages("devtools")
devtools::install_github("mothur/rchime")
```

## Usage

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
accepts [strollur](https://mothur.org/strollur/) objects or data.frames
as inputs. Let’s create a
[strollur::strollur](https://mothur.org/strollur/reference/strollur.html)
object using files from [mothur’s](https://mothur.org)
[Miseq_SOP](https://mothur.org/wiki/miseq_sop/) example analysis. Then
we will use the *denovo* method in
[`rchime()`](http://mothur.org/rchime/reference/rchime.md) to detect and
remove the chimeras from the dataset.

``` r

fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
abundance_data <- readRDS(rchime_example("miseq_abundance.rds"))

data <- strollur::new_dataset("rchime denovo example")

strollur::add(data, table = fasta_data, type = "sequences")
#> ℹ Added 6084 sequences.
#> [1] 6084
strollur::assign(data, table = abundance_data, type = "sequence_abundance")
#> ℹ Assigned 6084 sequence abundances.
#> [1] 6084

chimera_report <- rchime(data)
#> 
#> The denovo method runs with a single processor.
#> ℹ Added a chimera_report.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `4.26860308647156` seconds to detect and remove the chimeras.

data
#> rchime denovo example:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        3     0   2956.05
#> 25%-tile:        1  252    252      0        4     0  29551.50
#> Median:          1  252    252      0        4     0  59102.00
#> 75%-tile:        1  253    253      0        5     0  88652.50
#> 97.5%-tile:      1  253    253      0        6     0 115247.95
#> Maximum:         1  256    256      0        8     0 118202.00
#> Mean:            1  252    252      0        4     0      0.00
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime-chimeras   3588 10453
#> 
#> Number of unique seqs: 2496 
#> Total number of seqs: 118202 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1
```

## References

Many thanks for the great work of the
*[uchime](https://www.drive5.com/usearch/manual/chimeras.html)* and
*[vsearch](https://github.com/torognes/vsearch/graphs/contributors)*
teams!

- Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
  versatile open source tool for metagenomics. PeerJ 4:e2584. doi:
  10.7717/peerj.2584

- Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
  UCHIME improves sensitivity and speed of chimera detection.
  Bioinformatics 27:2194.
