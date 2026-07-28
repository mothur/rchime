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
#> 2.5%-tile:       1  252    252      0        4     0   3217.35
#> 25%-tile:        1  252    252      0        4     0  32164.50
#> Median:          1  253    253      0        4     0  64328.00
#> 75%-tile:        1  253    253      0        5     0  96491.50
#> 97.5%-tile:      1  254    254      0        6     0 125438.65
#> Maximum:         1  256    256      0        8     0 128655.00
#> Mean:            1  252    252      0        4     0  64328.00
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

rchime(strollur, reference = reference)
#> → rchime removed `1037` chimeras from your dataset.
#> → It took `12.0288896560669` seconds to detect and remove the chimeras.
strollur
#> rchime reference example:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   3191.43
#> 25%-tile:        1  252    252      0        4     0  31905.25
#> Median:          1  253    253      0        4     0  63809.50
#> 75%-tile:        1  253    253      0        5     0  95713.75
#> 97.5%-tile:      1  254    254      0        6     0 124427.57
#> Maximum:         1  256    256      0        8     0 127618.00
#> Mean:            1  252    252      0        4     0  63809.50
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime_chimeras    787  1037
#> 
#> Number of unique seqs: 5297 
#> Total number of seqs: 127618 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1

data <- rchime(data_df, reference = reference)
#> → rchime removed `1037` chimeras from your dataset.
#> → It took `12.0616800785065` seconds to detect and remove the chimeras.
data
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   3191.43
#> 25%-tile:        1  252    252      0        4     0  31905.25
#> Median:          1  253    253      0        4     0  63809.50
#> 75%-tile:        1  253    253      0        5     0  95713.75
#> 97.5%-tile:      1  254    254      0        6     0 124427.57
#> Maximum:         1  256    256      0        8     0 127618.00
#> Mean:            1  252    252      0        4     0  63809.50
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime_chimeras    787  1037
#> 
#> Number of unique seqs: 5297 
#> Total number of seqs: 127618 
#> 
#> Total number of custom reports: 1
```

## Results

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
returns a [strollur
object](https://mothur.org/strollur/reference/strollur.html). Let’s take
a closer look at the `chimera_report` generated by
[`rchime()`](http://mothur.org/rchime/reference/rchime.md). The
[chimera_report](https://mothur.org/rchime/articles/chimera_report.html)
is a data.frame with a row for each sequence in your dataset. Let’s take
a look at the first 5 chimeric sequences in the report:

``` r

chimera_report <- strollur::report(data, type = "chimera_report")

chimera_report[chimera_report$Chimeric_Status == "Y", ] |> head(n = 5)
#>        Score                                        Query          ParentA
#> 6  0.6724684  M00967_43_000000000-A3JHG_1_1101_10574_3283 7000004130085859
#> 9  0.4391696 M00967_43_000000000-A3JHG_1_1101_10956_16121       S000008023
#> 12 0.5777207 M00967_43_000000000-A3JHG_1_1101_11038_24430       S000008023
#> 20 1.0522959 M00967_43_000000000-A3JHG_1_1101_11538_27185       S000429249
#> 34 0.4824794 M00967_43_000000000-A3JHG_1_1101_13010_15146       S000149329
#>             ParentB       Top_Parent        QM       QA       QB      QAB
#> 6  7000004128190422 7000004130085859  93.22709 87.64940 81.67331 77.29084
#> 9  7000004131498097       S000008023  96.01594 93.22709 87.25100 83.66534
#> 12 7000004131498332       S000008023  96.41434 94.82072 80.87649 81.27490
#> 20 7000004131252252 7000004131252252 100.00000 82.47012 98.80478 81.27490
#> 34 7000004131498097       S000149329  91.23506 86.05578 78.08765 75.69721
#>          QT LY LN LA RY RN RA      Div Chimeric_Status
#> 6  87.64940 34  5 12 14  0  0 5.577689               Y
#> 9  93.22709 22  0  0 11  4  6 2.788845               Y
#> 12 94.82072 39  0  1  6  2  6 1.593625               Y
#> 20 98.80478  3  0  0 44  0  0 1.195219               Y
#> 34 86.05578 34  1 10 15  2  9 5.179283               Y
```
