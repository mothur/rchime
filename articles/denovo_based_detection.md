# Denovo Chimera Detection

## Overview

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
allows you to detect and remove chimeric sequences using the dataset as
it’s own reference (denovo). The denovo approach is our preferred method
for removing chimeras.
[`rchime()`](http://mothur.org/rchime/reference/rchime.md) can be used
with [strollur
objects](https://mothur.org/strollur/reference/strollur.html) or
data.frames as inputs. Let’s look at examples of both data types using
sequence data from [mothur’s](https://mothur.org)
[MiSeq_SOP](https://mothur.org/wiki/miseq_sop/) example analysis.

### Creating *[strollur](https://mothur.org/strollur)* objects

``` r
fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
abundance_data <- readRDS(rchime_example("miseq_abundance.rds"))

strollur <- strollur::new_dataset("rchime denovo example")

strollur::add(strollur, table = fasta_data, type = "sequences")
#> ℹ Added 6084 sequences.
#> [1] 6084
strollur::assign(strollur, table = abundance_data, type = "sequence_abundance")
#> ℹ Assigned 6084 sequence abundances.
#> [1] 6084

strollur
#> rchime denovo example:
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
df <- readRDS(rchime_example("miseq_data_frame_by_sample.rds"))

str(df)
#> 'data.frame':    11039 obs. of  4 variables:
#>  $ sequence_names: chr  "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10133_8460" ...
#>  $ sequences     : chr  "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGGAGTGCTTTTGAAACT"| __truncated__ ...
#>  $ samples       : chr  "F3D2" "F3D146" "F3D149" "F3D150" ...
#>  $ abundances    : int  222 1 1 1 127 17 32 13 95 86 ...
```

## Removing chimeras

When removing chimeras using the denovo method, the potential parents
are chosen from more abundant sequences in your dataset.

Before we remove the chimeras let’s discuss the `dereplicate` parameter.
If your dataset includes sample assignments, we recommend setting the
`dereplicate` parameter to `TRUE`. When `dereplicate` = `FALSE`, if a
sequence is flagged as chimeric in one sample, it is removed from all
samples. Our experience suggests that this is a bit aggressive since
we’ve seen rare sequences get flagged as chimeric when they’re the most
abundant sequence in another sample. For a more conservative approach,
set `dereplicate` = `TRUE` which will only remove sequences from the
samples in which they are flagged as chimeric. Let’s use the denovo
method to remove the chimeras.

``` r
strollur_results <- rchime(strollur, dereplicate = TRUE)
#> ℹ The denovo method runs with a single processor.
#> ℹ Added a chimera_report.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `7.60661172866821` seconds to detect and remove the chimeras.

strollur
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

data_frame_results <- rchime(df, dereplicate = TRUE)
#> ℹ The denovo method runs with a single processor.
#> → rchime detected `10453` chimeras in your dataset.
#> → It took `7.54418039321899` seconds to detect the chimeras.
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

### Chimeras

Results also contains a list of the names of the chimeric sequences.
Let’s get the names of the first 10 chimeras.

``` r
strollur_results$chimeras |> head(n = 10)
#>  [1] "M00967_43_000000000-A3JHG_1_1106_11629_14238"
#>  [2] "M00967_43_000000000-A3JHG_1_1103_26580_14708"
#>  [3] "M00967_43_000000000-A3JHG_1_2101_21700_24164"
#>  [4] "M00967_43_000000000-A3JHG_1_1112_23980_19089"
#>  [5] "M00967_43_000000000-A3JHG_1_2110_20944_24019"
#>  [6] "M00967_43_000000000-A3JHG_1_2107_23359_13368"
#>  [7] "M00967_43_000000000-A3JHG_1_1101_15516_19920"
#>  [8] "M00967_43_000000000-A3JHG_1_1103_19656_27166"
#>  [9] "M00967_43_000000000-A3JHG_1_1112_15051_19857"
#> [10] "M00967_43_000000000-A3JHG_1_2108_16511_16075"
```

### Set_abundance_values

Results will only contain the set_abundance_values list when dereplicate
= TRUE and you are not removing the chimeras automatically.
set_abundance_values has three items: ‘sequence_names’, ‘abundances’ and
‘samples’. For each sequence in your dataset there will be an entry in
sequence_names and abundances. The abundance values are parsed by
sample, and the order is given in set_abundance_values\$samples. Let’s
look at the first two sequences abundances after detecting the chimeras:

``` r
names(data_frame_results$set_abundance_values)
#> [1] "sequence_names" "abundances"     "samples"

data_frame_results$set_abundance_values$samples
#>  [1] "F3D0"   "F3D1"   "F3D141" "F3D142" "F3D143" "F3D144" "F3D145" "F3D146"
#>  [9] "F3D147" "F3D148" "F3D149" "F3D150" "F3D2"   "F3D3"   "F3D5"   "F3D6"  
#> [17] "F3D7"   "F3D8"   "F3D9"   "Mock"

sequences_names <- c(
  "M00967_43_000000000-A3JHG_1_1103_5171_14027",
  "M00967_43_000000000-A3JHG_1_1101_10133_8460"
)

df[df$sequence_names %in% sequences_names, c(1, 3, 4)]
#>                                   sequence_names samples abundances
#> 1    M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D2        222
#> 2    M00967_43_000000000-A3JHG_1_1101_10133_8460  F3D146          1
#> 3    M00967_43_000000000-A3JHG_1_1101_10133_8460  F3D149          1
#> 4    M00967_43_000000000-A3JHG_1_1101_10133_8460  F3D150          1
#> 5    M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D1        127
#> 6    M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D7         17
#> 7    M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D0         32
#> 8    M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D5         13
#> 9    M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D8         95
#> 10   M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D9         86
#> 11   M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D3          5
#> 12   M00967_43_000000000-A3JHG_1_1101_10133_8460    F3D6         20
#> 1092 M00967_43_000000000-A3JHG_1_1103_5171_14027    F3D0          5
#> 1093 M00967_43_000000000-A3JHG_1_1103_5171_14027  F3D148          1
#> 1094 M00967_43_000000000-A3JHG_1_1103_5171_14027    F3D2          4
#> 1095 M00967_43_000000000-A3JHG_1_1103_5171_14027    F3D1          6
#> 1096 M00967_43_000000000-A3JHG_1_1103_5171_14027  F3D142          1

data_frame_results$set_abundance_values$sequence_names[1:2]
#> [1] "M00967_43_000000000-A3JHG_1_1101_10133_8460" 
#> [2] "M00967_43_000000000-A3JHG_1_1101_10134_24617"

data_frame_results$set_abundance_values$abundances[1:2]
#> [[1]]
#>  [1]  32 127   0   0   0   0   0   1   0   0   1   1 222   5  13  20  17  95  86
#> [20]   0
#> 
#> [[2]]
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

We can see that the abundances for sequence
M00967_43_000000000-A3JHG_1_1101_10133_8460 remain the same meaning it
was not chimeric in any sample it was present in. We can see that
sequence M00967_43_000000000-A3JHG_1_1101_10134_24617 was found to be
chimeric in every sample it was included in. Now let’s look at an
example of a sequence that was found to be chimeric in some of the
samples it is present in.

``` r
df[
  df$sequence_names %in% "M00967_43_000000000-A3JHG_1_1103_5171_14027",
  c(1, 3, 4)
]
#>                                   sequence_names samples abundances
#> 1092 M00967_43_000000000-A3JHG_1_1103_5171_14027    F3D0          5
#> 1093 M00967_43_000000000-A3JHG_1_1103_5171_14027  F3D148          1
#> 1094 M00967_43_000000000-A3JHG_1_1103_5171_14027    F3D2          4
#> 1095 M00967_43_000000000-A3JHG_1_1103_5171_14027    F3D1          6
#> 1096 M00967_43_000000000-A3JHG_1_1103_5171_14027  F3D142          1

data_frame_results$set_abundance_values$samples
#>  [1] "F3D0"   "F3D1"   "F3D141" "F3D142" "F3D143" "F3D144" "F3D145" "F3D146"
#>  [9] "F3D147" "F3D148" "F3D149" "F3D150" "F3D2"   "F3D3"   "F3D5"   "F3D6"  
#> [17] "F3D7"   "F3D8"   "F3D9"   "Mock"

data_frame_results$set_abundance_values$abundances[
  data_frame_results$set_abundance_values$sequence_names %in%
    "M00967_43_000000000-A3JHG_1_1103_5171_14027"
]
#> [[1]]
#>  [1] 5 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

We can see that samples F3D0 and F3D142 did not find the sequence to be
chimeric, but F3D1, F3D2, F3D148 did find it to be chimeric so the
abundance for those samples is set to 0.
