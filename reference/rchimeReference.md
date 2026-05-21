# Detects chimeras using a reference based approach

detects chimeras in your data using a reference based approach.

## Usage

``` r
rchimeReference(
  sequence_name,
  sequence,
  reference_name,
  reference_sequence,
  options = NULL
)
```

## Arguments

- sequence_name, :

  vector of strings containing sequence names

- sequence, :

  vector of strings containing sequence nucleotide data

- reference_name, :

  vector of strings containing reference sequences names

- reference_sequence, :

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

results <- rchimeReference(sequence_name = fasta_data$sequence_name,
                           sequence = fasta_data$sequence,
                           reference_name = reference_data$sequence_name,
                           reference_sequence = reference_data$sequence)
```
