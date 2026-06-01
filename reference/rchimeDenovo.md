# Detects chimeras using a denovo method processing by sample.

detects chimeras from your data using a denovo method processing by
sample.

## Usage

``` r
rchimeDenovo(sequence_name, sequence, abundance, options = NULL)
```

## Arguments

- sequence_name, :

  2D vector of strings containing sequence names parsed by sample

- sequence, :

  2D vector of strings containing sequence nucleotide data parsed by
  sample

- abundance, :

  2D vector of sequence abundances parsed by sample

- options:

  list containing parameter options. Default = NULL.

## Value

list()

## See also

[`rchime()`](http://mothur.org/rchime/reference/rchime.md)

## Examples

``` r

# Read in data.frames containing sequence and abundance data

sequence_name <- readRDS(rchime_example("miseq_names_by_sample.rds"))
sequence <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
abundance <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))

# Detect chimeras

results <- rchimeDenovo(sequence_name = sequence_name,
                        sequence = sequence,
                        abundance = abundance)
```
