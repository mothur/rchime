# Detects chimeras using a denovo approach

detects chimeras in your data using a denovo approach.

## Usage

``` r
rchimeDenovoSingleSample(sequence_name, sequence, abundance, options = NULL)
```

## Arguments

- sequence_name, :

  vector of strings containing sequence names

- sequence, :

  vector of strings containing sequence nucleotide data

- abundance, :

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

results <- rchimeDenovoSingleSample(sequence_name = fasta_data$sequence_name,
                                    sequence = fasta_data$sequence,
                                    abundance = abundance_data$abundance)
```
