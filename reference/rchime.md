# Detect and remove chimeras from your [strollur](https://mothur.org/strollur/) dataset object or data.frame using a denovo approach or alternatively a reference model.

The `rchime()` function allows you to detect and remove chimeras from
your data using a denovo approach or alternatively a reference model.

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
