# Detect and remove chimeras from your [strollur](https://mothur.org/strollur/) dataset object.

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
allows you to detect and remove chimeras from your data using a denovo
approach or alternatively a reference model.

Our preferred way of doing this is to use the abundant sequences as our
reference (denovo).

This function uses code from the
[vsearch](https://github.com/torognes/vsearch) tools.

## Usage

``` r
# S3 method for class 'strollur'
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  silent = FALSE,
  remove_chimeras = TRUE,
  rchime_options = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence")
)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/) dataset object containing
  your sequence data.

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

- silent:

  Boolean, suppress console outputs. Default = FALSE.

- remove_chimeras:

  Boolean, remove chimeras from dataset. Default = TRUE.

- rchime_options:

  List, You can fine tune the vsearch specific options using the
  \[[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)\]
  function. Default = NULL.

- table_names, :

  Only used when `reference` is a data.frame. By default:

  table_names \<- list(sequence_name = "sequence_name", sequence =
  "sequence")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. Default column
  name is 'sequence_name'.

  In table_names, 'sequence' is a string containing the name of the
  column in 'table' that contains the sequences. Default column name is
  'sequence'.

## Value

list() containing a chimera report, and vector of the chimeric
sequence's names.

The [strollur](https://mothur.org/strollur/) dataset object will also be
updated. The sequences flagged as chimeric are removed and the chimera
report is added.

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

# Let's load a strollur object with 500 sequences

data <- strollur::load_dataset(
  rchime_example("strollur_multi_sample_small.rds")
)

# Detect and remove chimeras from the dataset using denovo approach by sample
# (recommended)

chimera_report <- rchime(data)
#> ℹ The denovo method runs with a single processor.
#> Added a chimera_report.
#> → rchime removed `128` chimeras from your dataset.
#> → It took `0.506592750549316` seconds to detect and remove the chimeras.
data
#> multi sample - 500 sequences:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  250    250      0        3     0    1.00
#> 2.5%-tile:       1  252    252      0        4     0  167.20
#> 25%-tile:        1  252    252      0        4     0 1672.00
#> Median:          1  253    253      0        4     0 3344.00
#> 75%-tile:        1  253    253      0        5     0 5016.00
#> 97.5%-tile:      1  253    253      0        6     0 6520.80
#> Maximum:         1  255    255      0        6     0 6688.00
#> Mean:            1  252    252      0        4     0 3344.14
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime-chimeras     90   128
#> 
#> Number of unique seqs: 410 
#> Total number of seqs: 6688 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1 
#> 
```
