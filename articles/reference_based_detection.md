# Reference Based Chimera Detection

## Overview

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
allows you to detect and remove chimeras from your dataset using a
referenced based approach.
[`rchime()`](http://mothur.org/rchime/reference/rchime.md) can be used
with [strollur
objects](https://mothur.org/strollur/reference/strollur.html) or
data.frames as inputs. Let’s look at examples of both data types using
sequence data from [mothur’s](https://mothur.org)
[MiSeq_SOP](https://mothur.org/wiki/miseq_sop/) example analysis, and
SILVA reference sequences trimmed to our region of interest as the
reference. *rchime* is designed to be flexible and you can use any
reference you choose.

### Creating *[strollur](https://mothur.org/strollur)* objects

Let’s create a [strollur
object](https://mothur.org/strollur/reference/strollur.html) that
contains our SILVA reference sequences.

``` r

reference <- strollur::new_dataset("Silva V4 Region")

strollur::add(reference,
  table = readRDS(rchime_example("reference.rds")),
  type = "sequence"
)
#> Added 14956 sequences.

reference
#> Silva V4 Region:
#> 
#>             starts ends nbases ambigs polymers numns  numseqs
#> Minimum:         1  270    270      0        3     0     1.00
#> 2.5%-tile:       1  292    292      0        4     0   374.90
#> 25%-tile:        1  293    293      0        4     0  3740.00
#> Median:          1  293    293      0        4     0  7479.00
#> 75%-tile:        1  293    293      0        5     0 11218.00
#> 97.5%-tile:      1  294    294      1        6     1 14583.10
#> Maximum:         1  351    351      5        9     5 14956.00
#> Mean:            1  292    292      0        4     0     0.00
#> 
#> Number of unique seqs: 14956 
#> Total number of seqs: 14956
```

Let’s create a [strollur
object](https://mothur.org/strollur/reference/strollur.html) containing
the MiSeq_SOP sequences.

``` r

fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
abundance_data <- readRDS(rchime_example("miseq_abundance.rds"))

strollur <- strollur::new_dataset("rchime reference example")

strollur::add(strollur, table = fasta_data, type = "sequence")
#> Added 6084 sequences.
strollur::assign(strollur, table = abundance_data, type = "sequence_abundance")
#> Assigned 6084 sequence abundances.

strollur
#> rchime reference example:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        3     0   3217.38
#> 25%-tile:        1  252    252      0        4     0  32164.75
#> Median:          1  252    252      0        4     0  64328.50
#> 75%-tile:        1  253    253      0        5     0  96492.25
#> 97.5%-tile:      1  253    253      0        6     0 125439.62
#> Maximum:         1  256    256      0        8     0 128655.00
#> Mean:            1  252    252      0        4     0      0.00
#> 
#> Number of unique seqs: 6084 
#> Total number of seqs: 128655 
#> 
#> Total number of samples: 20
```

### Loading data.frames

``` r

reference_data <- readRDS(rchime_example("reference.rds"))

str(reference_data)
#> 'data.frame':    14956 obs. of  2 variables:
#>  $ sequence_name: chr  "DQ415841.1" "L02888.1" "AB037561.1" "AF355050.1" ...
#>  $ sequence     : chr  "GTGCCAGCAGCCGCGGTAATACGGGGGGTGCAAGCGTTGTTCGGAATTATTGGGCGTAAAGAGCGTGTAGGCTGTCTGATGTGTCAGATGTGAAAGCCCTGGGCTTAACCC"| __truncated__ "GTGCCAGCAGCCGCGGTAATACGAGAGGGCTAGCGTTATTCGGAATTATTGGGCGTAAAGGGCGCGTAGGCTGGTTAATAAGTTAAAAGTGAAATCCCGAGGCTTAACCTT"| __truncated__ "GTGCCAGCAGCCGCGGTAATACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCT"| __truncated__ "GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATCACTGGGCGTAAAGGACGCGTAGGCGGGTTGATAAGTCAGATGTGAAATCCTACAGCTTAACTG"| __truncated__ ...

data_df <- readRDS(rchime_example("miseq_data_frame.rds"))

str(data_df)
#> 'data.frame':    6084 obs. of  3 variables:
#>  $ sequence_name: chr  "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10134_24617" "M00967_43_000000000-A3JHG_1_1101_10331_23332" "M00967_43_000000000-A3JHG_1_1101_10340_12294" ...
#>  $ sequence     : chr  "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGCAGGCGGCATGGCAAGTCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGCATTTGAAACT"| __truncated__ "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACGCCGTCGAGCCGTTGAAACTG"| __truncated__ ...
#>  $ abundance    : int  620 4 1 1 1 1 1 3 1 2 ...
```

## Removing chimeras

When removing chimeras using a reference, the potential parents are
chosen from the set of reference sequences. Let’s use the reference to
remove the chimeras.

``` r

strollur_results <- rchime(strollur, reference = reference)
#> Added a chimera_report.
#> → rchime removed `4849` chimeras from your dataset.
#> → It took `8.47708582878113` seconds to detect and remove the chimeras.

strollur
#> rchime reference example:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        3     0   3096.15
#> 25%-tile:        1  252    252      0        4     0  30952.50
#> Median:          1  252    252      0        4     0  61904.00
#> 75%-tile:        1  253    253      0        5     0  92855.50
#> 97.5%-tile:      1  253    253      0        6     0 120711.85
#> Maximum:         1  256    256      0        8     0 123806.00
#> Mean:            1  252    252      0        4     0      0.00
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence chimeras_rchime   1676  4849
#> 
#> Number of unique seqs: 4408 
#> Total number of seqs: 123806 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1

data_frame_results <- rchime(data_df, reference = reference_data)
#> → rchime detected `4849` chimeras in your dataset.
#> → It took `8.50946164131165` seconds to detect the chimeras.
```

## Results

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
returns a list containing the results of the function. When you are
running the command with a strollur object, the chimera_report is added,
and chimeras are removed for you automatically. Let’s take a closer look
at the results returned.

### Chimera Report

The
[chimera_report](https://mothur.org/rchime/articles/chimera_report.html)
is a data.frame with a row for each sequence in your dataset. Let’s take
a look at the first 5 chimeric sequences in the report:

``` r

strollur_results$chimera_report[
  strollur_results$chimera_report$Chimeric_Status == "Y",
] |> head(n = 5)
#>        Score                                        Query    ParentA    ParentB
#> 4  0.3599476 M00967_43_000000000-A3JHG_1_2113_29036_16812 AJ404681.1 AJ400237.1
#> 11 1.0186117  M00967_43_000000000-A3JHG_1_2108_19952_9563 AJ400267.1 AJ400236.1
#> 14 1.7641129  M00967_43_000000000-A3JHG_1_2107_6538_14332 AJ400267.1 AY230766.1
#> 16 2.7742347 M00967_43_000000000-A3JHG_1_1102_22525_13142 AJ400236.1 EU006128.1
#> 18 1.5735035 M00967_43_000000000-A3JHG_1_2108_19687_26893 AB021164.1 AJ400267.1
#>    Top_Parent        QM       QA       QB      QAB       QT LY LN LA RY RN RA
#> 4  AJ404681.1  94.82072 92.03187 86.85259 85.25896 92.03187 22  2 11  7  0  0
#> 11 AJ400267.1  98.80000 95.20000 91.60000 89.20000 95.20000 18  0  3  9  0  0
#> 14 AJ400267.1  97.93388 92.97521 77.68595 71.48760 92.97521 49  0  0 15  3  2
#> 16 AJ400236.1 100.00000 95.00000 87.91667 82.91667 95.00000 29  0  0 12  0  0
#> 18 AB021164.1  98.38710 93.14516 85.48387 80.64516 93.14516 33  1  0 13  0  3
#>         Div Chimeric_Status
#> 4  2.788845               Y
#> 11 3.600000               Y
#> 14 4.958678               Y
#> 16 5.000000               Y
#> 18 5.241935               Y
```

### Chimeras

Results also contains a list of the names of the chimeric sequences.
Let’s get the names of the first 10 chimeras.

``` r

strollur_results$chimeras |> head(n = 10)
#>  [1] "M00967_43_000000000-A3JHG_1_2113_29036_16812"
#>  [2] "M00967_43_000000000-A3JHG_1_2108_19952_9563" 
#>  [3] "M00967_43_000000000-A3JHG_1_2107_6538_14332" 
#>  [4] "M00967_43_000000000-A3JHG_1_1102_22525_13142"
#>  [5] "M00967_43_000000000-A3JHG_1_2108_19687_26893"
#>  [6] "M00967_43_000000000-A3JHG_1_2113_21308_17042"
#>  [7] "M00967_43_000000000-A3JHG_1_2101_14896_3554" 
#>  [8] "M00967_43_000000000-A3JHG_1_1109_10514_23344"
#>  [9] "M00967_43_000000000-A3JHG_1_1108_17549_23127"
#> [10] "M00967_43_000000000-A3JHG_1_1102_4170_17056"
```
