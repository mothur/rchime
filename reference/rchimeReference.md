# Detects chimeras using a reference based approach

detects chimeras in your data using a reference based approach.

## Usage

``` r
rchimeReference(
  sequence_names,
  sequences,
  reference_names,
  reference_sequences,
  options = NULL
)
```

## Arguments

- sequence_names, :

  vector of strings containing sequence names

- sequences, :

  vector of strings containing sequence nucleotide data

- reference_names, :

  vector of strings containing reference sequences names

- reference_sequences, :

  vector of strings containing reference sequences

- options:

  list containing parameter options. Default = NULL.

## See also

[`rchime()`](http://mothur.org/rchime/reference/rchime.md)

## Examples

``` r
# Read in data.frames containing sequence data

fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
reference_data <- readRDS(rchime_example("reference.rds"))

# Detect chimeras

results <- rchimeReference(sequence_names = fasta_data$sequence_names,
                           sequences = fasta_data$sequences,
                           reference_names = reference_data$sequence_names,
                           reference_sequences = reference_data$sequences)
```
