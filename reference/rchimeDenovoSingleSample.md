# Detects chimeras using a denovo approach

detects chimeras in your data using a denovo approach.

## Usage

``` r
rchimeDenovoSingleSample(sequence_names, sequences, abundances, options = NULL)
```

## Arguments

- sequence_names, :

  vector of strings containing sequence names

- sequences, :

  vector of strings containing sequence nucleotide data

- abundances, :

  vector of sequence abundances

- options:

  list containing parameter options. Default = NULL.

## See also

[`rchime()`](http://mothur.org/rchime/reference/rchime.md)

## Examples

``` r
# Read in data.frames containing sequence and abundance data

fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
abundance_data <- readRDS(rchime_example("single_sample_abundance.rds"))

# Detect chimeras

results <- rchimeDenovoSingleSample(sequence_names = fasta_data$sequence_names,
                                    sequences = fasta_data$sequences,
                                    abundances = abundance_data$abundances)
```
