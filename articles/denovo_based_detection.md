# De novo Chimera Detection

## Overview

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
allows you to detect and remove chimeric sequences using the dataset as
it’s own reference (de novo). The denovo approach is our preferred
method for removing chimeras.
[`rchime()`](http://mothur.org/rchime/reference/rchime.md) can be used
with [strollur
objects](https://mothur.org/strollur/reference/strollur.html) or
data.frames as inputs. Let’s look at examples of both data types using
sequence data from [mothur’s](https://mothur.org)
[MiSeq_SOP](https://mothur.org/wiki/miseq_sop/) example analysis.

### Creating *[strollur](https://mothur.org/strollur/)* objects

``` r

fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
abundance_data <- readRDS(rchime_example("miseq_abundance.rds"))

strollur <- strollur::new_dataset("rchime de novo example")

strollur::add(strollur, table = fasta_data, type = "sequence")
#> Added 6084 sequences.
strollur::assign(strollur, table = abundance_data, type = "sequence_abundance")
#> Assigned 6084 sequence abundances.

strollur
#> rchime de novo example:
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

df <- readRDS(rchime_example("miseq_data_frame_by_sample.rds"))

str(df)
#> 'data.frame':    11039 obs. of  4 variables:
#>  $ sequence_name: chr  "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10133_8460" ...
#>  $ sequence     : chr  "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ ...
#>  $ sample       : chr  "F3D2" "F3D146" "F3D149" "F3D150" ...
#>  $ abundance    : int  222 1 1 1 127 17 32 13 95 86 ...
```

## Removing chimeras

When removing chimeras using the de novo method, the potential parents
are chosen from more abundant sequences in your dataset.

Before we remove the chimeras let’s discuss the `dereplicate` parameter.
When `dereplicate=FALSE`, if a sequence is flagged as chimeric in one
sample, it is removed from all samples. Our experience suggests that
this is a bit aggressive since we’ve seen rare sequences get flagged as
chimeric when they’re the most abundant sequence in another sample. For
a more conservative approach, we recommend using the default
`dereplicate=TRUE` which will only remove sequences from the samples in
which they are flagged as chimeric. Let’s use the de novo method to
remove the chimeras.

``` r

rchime(strollur)
#> ℹ The de novo method runs with a single processor.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `7.50702977180481` seconds to detect and remove the chimeras.
strollur
#> rchime de novo example:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   2956.03
#> 25%-tile:        1  252    252      0        4     0  29551.25
#> Median:          1  253    253      0        4     0  59101.50
#> 75%-tile:        1  253    253      0        5     0  88651.75
#> 97.5%-tile:      1  254    254      0        6     0 115246.97
#> Maximum:         1  256    256      0        8     0 118202.00
#> Mean:            1  252    252      0        4     0  59101.50
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime_chimeras   3588 10453
#> 
#> Number of unique seqs: 2496 
#> Total number of seqs: 118202 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1

data <- rchime(df)
#> ℹ The de novo method runs with a single processor.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `7.48063826560974` seconds to detect and remove the chimeras.
data
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   2956.03
#> 25%-tile:        1  252    252      0        4     0  29551.25
#> Median:          1  253    253      0        4     0  59101.50
#> 75%-tile:        1  253    253      0        5     0  88651.75
#> 97.5%-tile:      1  254    254      0        6     0 115246.97
#> Maximum:         1  256    256      0        8     0 118202.00
#> Mean:            1  252    252      0        4     0  59101.50
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime_chimeras   3588 10453
#> 
#> Number of unique seqs: 2496 
#> Total number of seqs: 118202 
#> 
#> Total number of samples: 20 
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
#>        Score                                        Query
#> 66 0.3805621 M00967_43_000000000-A3JHG_1_1106_11629_14238
#> 70 0.5261480 M00967_43_000000000-A3JHG_1_1103_26580_14708
#> 80 0.5580357 M00967_43_000000000-A3JHG_1_2101_21700_24164
#> 89 0.3348214 M00967_43_000000000-A3JHG_1_1112_23980_19089
#> 91 0.6377551 M00967_43_000000000-A3JHG_1_2110_20944_24019
#>                                         ParentA
#> 66 M00967_43_000000000-A3JHG_1_1107_15750_18592
#> 70 M00967_43_000000000-A3JHG_1_2110_12856_16229
#> 80 M00967_43_000000000-A3JHG_1_1113_11294_24024
#> 89 M00967_43_000000000-A3JHG_1_1107_15750_18592
#> 91 M00967_43_000000000-A3JHG_1_2101_22400_13416
#>                                         ParentB
#> 66 M00967_43_000000000-A3JHG_1_2101_22400_13416
#> 70 M00967_43_000000000-A3JHG_1_1107_15750_18592
#> 80 M00967_43_000000000-A3JHG_1_2110_12856_16229
#> 89 M00967_43_000000000-A3JHG_1_1109_25348_18015
#> 91  M00967_43_000000000-A3JHG_1_1112_6862_18037
#>                                      Top_Parent        QM       QA       QB
#> 66 M00967_43_000000000-A3JHG_1_1107_15750_18592  99.60317 98.01587 94.44444
#> 70 M00967_43_000000000-A3JHG_1_2110_12856_16229 100.00000 97.61905 95.63492
#> 80 M00967_43_000000000-A3JHG_1_1113_11294_24024 100.00000 98.01587 94.44444
#> 89 M00967_43_000000000-A3JHG_1_1107_15750_18592 100.00000 98.80952 94.44444
#> 91 M00967_43_000000000-A3JHG_1_2101_22400_13416 100.00000 98.01587 93.65079
#>         QAB       QT LY LN LA RY RN RA      Div Chimeric_Status
#> 66 93.25397 98.01587 13  0  0  4  0  1 1.587302               Y
#> 70 93.25397 97.61905 11  0  0  6  0  0 2.380952               Y
#> 80 92.46032 98.01587 14  0  0  5  0  0 1.984127               Y
#> 89 93.25397 98.80952 14  0  0  3  0  0 1.190476               Y
#> 91 91.66667 98.01587 16  0  0  5  0  0 1.984127               Y
```
