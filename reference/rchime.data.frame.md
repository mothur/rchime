# Detect chimeras in your data.frames.

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
allows you to detect chimeras from your data using a denovo approach or
alternatively a reference model.

Our preferred way of doing this is to use the abundant sequences as our
reference (denovo).

This function uses code from the
[vsearch](https://github.com/torognes/vsearch) tools.

## Usage

``` r
# S3 method for class 'data.frame'
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  verbose = TRUE,
  remove_chimeras = NULL,
  rchime_options = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence", abundance =
    "abundance", sample = "sample")
)
```

## Arguments

- data:

  a data.frame containing your sequence data.

- reference:

  a data.frame containing reference sequence data.

- dereplicate:

  logical. The dereplicate option allows you to flag chimeras by sample.
  When `dereplicate=FALSE`, if a sequence is flagged as chimeric in one
  sample, it is removed from all samples. Our experience suggests that
  this is a bit aggressive since we’ve seen rare sequences get flagged
  as chimeric when they’re the most abundant sequence in another sample.
  For a more conservative approach, we recommend using the default
  `dereplicate=TRUE` which will only remove sequences from the samples
  in which they are flagged as chimeric.

- verbose:

  logical, allow console outputs. Default = `TRUE`.

- remove_chimeras:

  Only used when `data` is a strollur object.

- rchime_options:

  List, You can fine tune the vsearch specific options using the
  \[[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)\]
  function. Default = NULL.

- table_names, :

  named list used to indicate the names of the columns in the
  data.frame. Only used when `data` is a data.frame. By default:

  table_names \<- list(sequence_name = "sequence_name", sequence =
  "sequence" abundance = "abundance", sample = "sample")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. Default column
  name is 'sequence_name'.

  In table_names, 'sequence' is a string containing the name of the
  column in 'table' that contains the sequences. Default column name is
  'sequence'.

  In table_names, 'abundance' is a string containing the name of the
  column in 'table' that contains the abundances. Default column name is
  'abundance'.

  In table_names, 'sample' is a string containing the name of the column
  in 'table' that contains the samples. Default column name is 'sample'.

## Value

list() containing a chimera report, and vector of the chimeric
sequence's names. If running with multiple samples and dereplicate =
TRUE, then a table containing the modified sequence abundances will also
be returned.

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

# Detect chimeras from the dataset using denovo approach by sample
# (recommended)

data <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))

chimera_report <- rchime(data)
#> ℹ The denovo method runs with a single processor.
#> → rchime detected `128` chimeras in your dataset.
#> → It took `0.543200969696045` seconds to detect the chimeras.
```
