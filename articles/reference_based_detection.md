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
the [`silva_gold()`](http://mothur.org/rchime/reference/silva_gold.md)
reference sequences. *rchime* is designed to be flexible and you can use
any reference you choose.

### Creating *[strollur](https://mothur.org/strollur/)* objects

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
#> 2.5%-tile:       1  252    252      0        4     0   3216.38
#> 25%-tile:        1  252    252      0        4     0  32163.75
#> Median:          1  253    253      0        4     0  64327.50
#> 75%-tile:        1  253    253      0        5     0  96491.25
#> 97.5%-tile:      1  254    254      0        6     0 125438.62
#> Maximum:         1  256    256      0        8     0 128655.00
#> Mean:            1  252    252      0        4     0  64327.64
#> 
#> Number of unique seqs: 6084 
#> Total number of seqs: 128655 
#> 
#> Total number of samples: 20
```

### Loading data.frames

``` r

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

reference <- silva_gold()
str(reference)
#> 'data.frame':    5181 obs. of  2 variables:
#>  $ sequence_name: chr  "7000004128189528" "7000004128189537" "7000004128189547" "7000004128189554" ...
#>  $ sequence     : chr  "GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGAAAGGCCCTTCGGGGTACTCGAGCGGCGAACGGGTGAGTAACACGTGGGCAACCTACCCCCAGCACCGG"| __truncated__ "GATGAACGCTGGCGGTATGCTTAACACATGCAAGTCGAACGGAATCTTCGGATTTAGTGGCGGACGGGTGAGTAACGCGTGAGAATCTAGCTCTAGGTCGGGGACAACCAC"| __truncated__ "ATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGTAACAGGTCTTCGGATGCTGACGAGTGGCGAACGGGTGAGTAATACATCGGAACGTGCCCGATCGTGGG"| __truncated__ "GATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGAAGCATCTTCGGATGCTTAGTGGCGAACGGGTGAGTAACACGTAGATAACCTACCTTTAACTCGAGGATA"| __truncated__ ...

strollur_results <- rchime(strollur, reference = reference)
#> Added a chimera_report report.
#> → rchime removed `1037` chimeras from your dataset.
#> → It took `11.683434009552` seconds to detect and remove the chimeras.

strollur
#> rchime reference example:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   3190.45
#> 25%-tile:        1  252    252      0        4     0  31904.50
#> Median:          1  253    253      0        4     0  63809.00
#> 75%-tile:        1  253    253      0        5     0  95713.50
#> 97.5%-tile:      1  254    254      0        6     0 124427.55
#> Maximum:         1  256    256      0        8     0 127618.00
#> Mean:            1  252    252      0        4     0  63809.14
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence chimeras_rchime    787  1037
#> 
#> Number of unique seqs: 5297 
#> Total number of seqs: 127618 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1

data_frame_results <- rchime(data_df, reference = reference)
#> → rchime detected `1037` chimeras in your dataset.
#> → It took `11.727872133255` seconds to detect the chimeras.
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
#>        Score                                        Query          ParentA
#> 14 0.3524927  M00967_43_000000000-A3JHG_1_2107_6538_14332       S000013923
#> 18 0.4054292 M00967_43_000000000-A3JHG_1_2108_19687_26893       S000008023
#> 21 0.3314002 M00967_43_000000000-A3JHG_1_2113_21308_17042       S000015682
#> 23 0.3266551  M00967_43_000000000-A3JHG_1_2110_13110_1977       S000260075
#> 30 0.6744335 M00967_43_000000000-A3JHG_1_1109_10514_23344 7000004128490675
#>             ParentB       Top_Parent       QM       QA       QB      QAB
#> 14       S000022285       S000022285 89.03509 80.70175 82.01754 77.63158
#> 18 7000004128504291       S000008023 94.82072 93.22709 77.29084 77.68924
#> 21       S000017014       S000015682 97.98387 96.37097 85.88710 84.67742
#> 23 7000004128191216       S000260075 95.93496 94.30894 79.26829 76.42276
#> 30 7000004131498332 7000004128490675 94.42231 88.84462 77.68924 74.10359
#>          QT LY LN LA RY RN RA      Div Chimeric_Status
#> 14 82.01754 22  6 14 19  0  5 7.017544               Y
#> 18 93.22709 46  2 10  4  0  1 1.593625               Y
#> 21 96.37097 32  2  0  4  0  3 1.612903               Y
#> 23 94.30894 45  4  6  4  0  0 1.626016               Y
#> 30 88.84462 45  3  3 15  1  7 5.577689               Y
```

### Chimeras

Results also contains a list of the names of the chimeric sequences.
Let’s get the names of the first 10 chimeras.

``` r

strollur_results$chimeras |> head(n = 10)
#>  [1] "M00967_43_000000000-A3JHG_1_2107_6538_14332" 
#>  [2] "M00967_43_000000000-A3JHG_1_2108_19687_26893"
#>  [3] "M00967_43_000000000-A3JHG_1_2113_21308_17042"
#>  [4] "M00967_43_000000000-A3JHG_1_2110_13110_1977" 
#>  [5] "M00967_43_000000000-A3JHG_1_1109_10514_23344"
#>  [6] "M00967_43_000000000-A3JHG_1_2111_20856_13433"
#>  [7] "M00967_43_000000000-A3JHG_1_1110_15964_18711"
#>  [8] "M00967_43_000000000-A3JHG_1_2101_14761_24747"
#>  [9] "M00967_43_000000000-A3JHG_1_1114_20514_18980"
#> [10] "M00967_43_000000000-A3JHG_1_1113_26305_21260"
```
