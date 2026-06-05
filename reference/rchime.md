# Detect and remove chimeras from your [strollur](https://mothur.org/strollur/) object or data.frame using a denovo approach or alternatively a reference model.

The `rchime()` function allows you to detect and remove chimeras from
your data using a denovo approach or alternatively a reference model.

Our preferred way of doing this is to use the abundant sequences as our
reference (denovo).

This function uses code from the
[vsearch](https://github.com/torognes/vsearch) tools.

## Usage

``` r
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  silent = FALSE,
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

- silent:

  Boolean, suppress console outputs. Default = FALSE.

- remove_chimeras:

  Boolean, remove chimeras from dataset. Default = TRUE. Only used when
  `data` is a strollur object.

- rchime_options:

  List, You can fine tune the vsearch specific options using the
  \[[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)\]
  function. Default = NULL.

- table_names, :

  named list used to indicate the names of the columns in the
  data.frame. Only used when `data` is a data.frame.

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

# Let's use a strollur object with 100 sequences

data <- strollur::load_dataset(
  rchime_example("strollur_multi_sample_tiny.rds")
)

chimera_report <- rchime(data)
#> ℹ The denovo method runs with a single processor.
#> Added a chimera_report.
#> → rchime removed `1` chimeras from your dataset.
#> → It took `0.186153411865234` seconds to detect and remove the chimeras.
data
#> rchime denovo example - 100 sequences:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  251    251      0        3     0    1.00
#> 2.5%-tile:       1  251    251      0        4     0   51.68
#> 25%-tile:        1  252    252      0        4     0  516.75
#> Median:          1  253    253      0        5     0 1033.50
#> 75%-tile:        1  253    253      0        5     0 1550.25
#> 97.5%-tile:      1  253    253      0        6     0 2015.33
#> Maximum:         1  255    255      0        8     0 2067.00
#> Mean:            1  252    252      0        4     0 1033.64
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime-chimeras      1     1
#> 
#> Number of unique seqs: 99 
#> Total number of seqs: 2067 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1 
#> 
```
