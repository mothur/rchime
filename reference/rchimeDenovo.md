# Detects chimeras using a denovo method processing by sample.

detects chimeras from your data using a denovo method processing by
sample.

## Usage

``` r
rchimeDenovo(sequence_names, sequences, abundances, options = NULL)
```

## Arguments

- sequence_names, :

  2D vector of strings containing sequence names parsed by sample

- sequences, :

  2D vector of strings containing sequence nucleotide data parsed by
  sample

- abundances, :

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

sequence_names <- readRDS(rchime_example("miseq_names_by_sample.rds"))
sequences <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
abundances <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))

# Detect chimeras

options <- rchime_options(dereplicate = TRUE)
results <- rchimeDenovo(sequence_names = sequence_names,
                        sequences = sequences,
                        abundances = abundances,
                        options)
```
