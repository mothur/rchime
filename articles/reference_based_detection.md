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

strollur_results <- rchime(strollur, reference = reference)
#> → rchime removed `1037` chimeras from your dataset.
#> → It took `11.7914447784424` seconds to detect and remove the chimeras.

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

data_frame_results <- rchime(data_df, reference = reference)
#> → rchime removed `1037` chimeras from your dataset.
#> → It took `11.8895835876465` seconds to detect and remove the chimeras.
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
