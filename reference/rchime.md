# rchime

The `rchime()` function allows you to detect and remove chimeras from
your data using a de novo approach or alternatively a reference model.

Our preferred way of doing this is to use the abundant sequences as our
reference (de novo).

This function uses code from the
[vsearch](https://github.com/torognes/vsearch) tools.

## Usage

``` r
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  verbose = TRUE,
  remove_chimeras = TRUE,
  rchime_options = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence", abundance =
    "abundance", sample = "sample")
)

# S3 method for class 'strollur'
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  verbose = TRUE,
  remove_chimeras = TRUE,
  rchime_options = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence")
)

# S3 method for class 'data.frame'
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  verbose = TRUE,
  remove_chimeras = TRUE,
  rchime_options = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence", abundance =
    "abundance", sample = "sample")
)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/) dataset object or a
  data.frame containing your sequence data.

- reference:

  a [strollur](https://mothur.org/strollur/) dataset object or a
  data.frame containing reference sequence data.

- dereplicate:

  logical. The dereplicate option allows you to remove chimeras by
  sample. When `dereplicate=FALSE`, if a sequence is flagged as chimeric
  in one sample, it is removed from all samples. Our experience suggests
  that this is a bit aggressive since we’ve seen rare sequences get
  flagged as chimeric when they’re the most abundant sequence in another
  sample. For a more conservative approach, we recommend using the
  default `dereplicate=TRUE` which will only remove sequences from the
  samples in which they are flagged as chimeric.

- verbose:

  logical, allow console outputs. Default = `TRUE`.

- remove_chimeras:

  Boolean, remove chimeras from dataset. Default = `TRUE`.

- rchime_options:

  List, You can fine tune the vsearch specific options using the
  \[[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)\]
  function. Default = `NULL`.

- table_names, :

  named list used to indicate the names of the columns in the
  data.frame. Only used when `data` or `reference` are a data.frames.
  Default = table_names \<- list(sequence_name = "sequence_name",
  sequence = "sequence")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. Default column
  name is 'sequence_name'.

  In table_names, 'sequence' is a string containing the name of the
  column in 'table' that contains the sequences. Default column name is
  'sequence'.

## Value

data a [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## References

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
versatile open source tool for metagenomics. PeerJ 4:e2584. doi:
10.7717/peerj.2584

Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
UCHIME improves sensitivity and speed of chimera detection.
Bioinformatics 27:2194.

## See also

[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)
to set vsearch specific parameters.

## Author

Sarah Westcott, <swestcot@umich.edu>

## Examples

``` r

# Let's use a strollur object with 500 sequences

data <- strollur::load_dataset(
  rchime_example("strollur_multi_sample_small.rds")
)

rchime(data)
#> ℹ The de novo method runs with a single processor.
#> → rchime removed `128` chimeras from your dataset.
#> → It took `0.536090135574341` seconds to detect and remove the chimeras.
data
#> miseq - 500 sequences:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  250    250      0        3     0    1.00
#> 2.5%-tile:       1  252    252      0        4     0  168.18
#> 25%-tile:        1  252    252      0        4     0 1672.75
#> Median:          1  253    253      0        4     0 3344.50
#> 75%-tile:        1  253    253      0        5     0 5016.25
#> 97.5%-tile:      1  253    253      0        6     0 6520.82
#> Maximum:         1  255    255      0        6     0 6688.00
#> Mean:            1  252    252      0        4     0 3344.50
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime_chimeras     90   128
#> 
#> Number of unique seqs: 410 
#> Total number of seqs: 6688 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1 
#> 


# Detect and remove chimeras from the dataset using de novo approach by
# sample (recommended)

table <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))

data <- rchime(data = table)
#> ℹ The de novo method runs with a single processor.
#> → rchime removed `128` chimeras from your dataset.
#> → It took `0.523041963577271` seconds to detect and remove the chimeras.
data
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  250    250      0        3     0    1.00
#> 2.5%-tile:       1  252    252      0        4     0  168.18
#> 25%-tile:        1  252    252      0        4     0 1672.75
#> Median:          1  253    253      0        4     0 3344.50
#> 75%-tile:        1  253    253      0        5     0 5016.25
#> 97.5%-tile:      1  253    253      0        6     0 6520.82
#> Maximum:         1  255    255      0        6     0 6688.00
#> Mean:            1  252    252      0        4     0 3344.50
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime_chimeras     90   128
#> 
#> Number of unique seqs: 410 
#> Total number of seqs: 6688 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1 
#> 
```
