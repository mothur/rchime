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

strollur_results <- rchime(strollur)
#> ℹ The de novo method runs with a single processor.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `7.46920990943909` seconds to detect and remove the chimeras.

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

data_frame_results <- rchime(df)
#> ℹ The de novo method runs with a single processor.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `7.46071720123291` seconds to detect and remove the chimeras.
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
#> NULL
```

### Chimeras

Results also contains a list of the names of the chimeric sequences.
Let’s get the names of the first 10 chimeras.

``` r

strollur_results$chimeras |> head(n = 10)
#> NULL
```

### Set_abundance_values

Results will only contain the set_abundance_values list when dereplicate
= TRUE and you are NOT removing the chimeras automatically.
set_abundance_values has three items: ‘sequence_names’, ‘abundances’ and
‘samples’. For each sequence in your dataset there will be an entry in
sequence_names and abundances. The abundance values are parsed by
sample, and the order is given in set_abundance_values\$samples. Let’s
look at the first two sequences abundances after detecting the chimeras:

``` r

names(data_frame_results$set_abundance_values)
#> NULL

data_frame_results$set_abundance_values$samples
#> NULL

sequences_names <- c(
  "M00967_43_000000000-A3JHG_1_1103_5171_14027",
  "M00967_43_000000000-A3JHG_1_1101_10133_8460"
)

df[df$sequence_name %in% sequences_names, c(1, 3, 4)]
#>                                    sequence_name sample abundance
#> 1    M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D2       222
#> 2    M00967_43_000000000-A3JHG_1_1101_10133_8460 F3D146         1
#> 3    M00967_43_000000000-A3JHG_1_1101_10133_8460 F3D149         1
#> 4    M00967_43_000000000-A3JHG_1_1101_10133_8460 F3D150         1
#> 5    M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D1       127
#> 6    M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D7        17
#> 7    M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D0        32
#> 8    M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D5        13
#> 9    M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D8        95
#> 10   M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D9        86
#> 11   M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D3         5
#> 12   M00967_43_000000000-A3JHG_1_1101_10133_8460   F3D6        20
#> 1092 M00967_43_000000000-A3JHG_1_1103_5171_14027   F3D0         5
#> 1093 M00967_43_000000000-A3JHG_1_1103_5171_14027 F3D148         1
#> 1094 M00967_43_000000000-A3JHG_1_1103_5171_14027   F3D2         4
#> 1095 M00967_43_000000000-A3JHG_1_1103_5171_14027   F3D1         6
#> 1096 M00967_43_000000000-A3JHG_1_1103_5171_14027 F3D142         1

data_frame_results$set_abundance_values$sequence_name[1:2]
#> NULL

data_frame_results$set_abundance_values$abundance[1:2]
#> NULL
```

We can see that the abundances for sequence
M00967_43_000000000-A3JHG_1_1101_10133_8460 remain the same meaning it
was not chimeric in any sample it was present in. We can see that
sequence M00967_43_000000000-A3JHG_1_1101_10134_24617 was found to be
chimeric in every sample it was included in.

Now let’s look at an example of a sequence that was found to be chimeric
in some of the samples it is present in.

``` r

df[
  df$sequence_name %in% "M00967_43_000000000-A3JHG_1_1103_5171_14027",
  c(1, 3, 4)
]
#>                                    sequence_name sample abundance
#> 1092 M00967_43_000000000-A3JHG_1_1103_5171_14027   F3D0         5
#> 1093 M00967_43_000000000-A3JHG_1_1103_5171_14027 F3D148         1
#> 1094 M00967_43_000000000-A3JHG_1_1103_5171_14027   F3D2         4
#> 1095 M00967_43_000000000-A3JHG_1_1103_5171_14027   F3D1         6
#> 1096 M00967_43_000000000-A3JHG_1_1103_5171_14027 F3D142         1

data_frame_results$set_abundance_values$samples
#> NULL

data_frame_results$set_abundance_values$abundance[
  data_frame_results$set_abundance_values$sequence_name %in%
    "M00967_43_000000000-A3JHG_1_1103_5171_14027"
]
#> NULL
```

We can see that samples F3D0 and F3D142 did not find the sequence to be
chimeric, but F3D1, F3D2, F3D148 did find it to be chimeric so the
abundance for those samples is set to 0.
