# Detect and remove chimeras from your [strollur](https://mothur.org/strollur/) dataset object.

The [`rchime()`](http://mothur.org/rchime/reference/rchime.md) function
allows you to detect and remove chimeras from your data using a denovo
approach or alternatively a reference model.

Our preferred way of doing this is to use the abundant sequences as our
reference (denovo). In addition, if a sequence is flagged as chimeric in
one sample, and (dereplicate=FALSE) then sequence is removed from all
samples. Our experience suggests that this is a bit aggressive since
we’ve seen rare sequences get flagged as chimeric when they’re the most
abundant sequence in another sample. For a more conservative approach,
set (dereplicate=TRUE) which will only remove sequences from the samples
in which they are flagged as chimeric.

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
  table_names = NULL
)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/) dataset object containing
  your sequence data.

- reference:

  a [strollur](https://mothur.org/strollur/) dataset object containing
  reference sequence data.

- dereplicate, :

  Boolean, The dereplicate option allows you to remove chimeras by
  sample. For example, if dereplicate parameter is FALSE, then if one
  group finds the sequence to be chimeric, it will be removed from all
  groups. If dereplicate is set to TRUE, sequences found to be chimeric
  are only removed from the sample they are found to be chimeric in.
  Default = TRUE.

- silent:

  Boolean, suppress console outputs. Default = FALSE.

- remove_chimeras:

  Boolean, remove chimeras from dataset. Default = TRUE.

- rchime_options:

  List, You can fine tune the vsearch specific options using the
  \[[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)\]
  function. Default = NULL.

- table_names, :

  Only used when `data` is a data.frame. Default = NULL.

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

# Let's load a strollur object named "rchime denovo example"

data_denovo <- strollur::load_dataset(
  rchime_example("strollur_multi_sample.rds")
)

# Detect and remove chimeras from the dataset using denovo approach by sample
# (recommended)

chimera_report <- rchime(data_denovo, dereplicate = TRUE)
#> ℹ The denovo method runs with a single processor.
#> Added a chimera_report.
#> → rchime removed `10453` chimeras from your dataset.
#> → It took `7.76756072044373` seconds to detect and remove the chimeras.
data_denovo
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   2956.05
#> 25%-tile:        1  252    252      0        4     0  29551.50
#> Median:          1  253    253      0        4     0  59102.00
#> 75%-tile:        1  253    253      0        5     0  88652.50
#> 97.5%-tile:      1  254    254      0        6     0 115247.95
#> Maximum:         1  256    256      0        8     0 118203.00
#> Mean:            1  252    252      0        4     0  59102.00
#> 
#> scrap_summary:
#>       type      trash_code unique total
#> 1 sequence rchime-chimeras   3588 10453
#> 
#> Number of unique seqs: 2496 
#> Total number of seqs: 118202 
#> 
#> Total number of samples: 20 
#> Total number of custom reports: 1 
#> 

# Alternatively you can detect and remove chimeras using a reference

data_reference <- strollur::load_dataset(
  rchime_example("strollur_multi_sample.rds")
)

reference <- strollur::load_dataset(rchime_example("strollur_reference.rds"))
#> Error: std::bad_alloc

chimera_report <- rchime(data_reference, reference = reference)
#> Error: object 'reference' not found
data_reference
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  249    249      0        3     0      1.00
#> 2.5%-tile:       1  252    252      0        4     0   3217.38
#> 25%-tile:        1  252    252      0        4     0  32164.75
#> Median:          1  253    253      0        4     0  64328.50
#> 75%-tile:        1  253    253      0        5     0  96492.25
#> 97.5%-tile:      1  254    254      0        6     0 125439.62
#> Maximum:         1  256    256      0        8     0 128656.00
#> Mean:            1  252    252      0        4     0  64328.50
#> 
#> Number of unique seqs: 6084 
#> Total number of seqs: 128655 
#> 
#> Total number of samples: 20 
#> 
```
