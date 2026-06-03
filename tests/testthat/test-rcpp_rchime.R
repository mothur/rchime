# test rpp_rchime

test_that("test rchimeReference - errors", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  # not equal sizes
  expect_error(rchimeReference(
    sequence_name = fasta_data$sequence_name[1:50],
    sequence = fasta_data$sequence[1:100],
    reference_name = reference_data$sequence_name,
    reference_sequence = reference_data$sequence
  ))

  # not equal sizes
  expect_error(rchimeReference(
    sequence_name = fasta_data$sequence_name[1:100],
    sequence = fasta_data$sequence[1:100],
    reference_name = reference_data$sequence_name[1:5],
    reference_sequence = reference_data$sequence[1:10]
  ))

  options <- list(xn = 8.0)
  # create options but leave off processors and dereplicate
  results <- rchimeReference(
    sequence_name = fasta_data$sequence_name[1:100],
    sequence = fasta_data$sequence[1:100],
    reference_name = reference_data$sequence_name,
    reference_sequence = reference_data$sequence,
    options = options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 24)
  expect_equal(nrow(results$chimera_report), 100)
})

test_that("test rchimeDenovo, rchimeDenovoSingleSample - errors", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  abundance_data <- readRDS(rchime_example("single_sample_abundance.rds"))

  # not equal sizes
  expect_error(rchimeDenovoSingleSample(
    sequence_name = fasta_data$sequence_name[1:50],
    sequence = fasta_data$sequence[1:100],
    abundance = abundance_data$abundance
  ))

  # not equal sizes
  expect_error(rchimeDenovoSingleSample(
    sequence_name = fasta_data$sequence_name[1:100],
    sequence = fasta_data$sequence[1:100],
    abundance = abundance_data$abundance
  ))

  options <- list(xn = 8.0)
  # create options but leave off processors and dereplicate
  results <- rchimeDenovoSingleSample(
    sequence_name = fasta_data$sequence_name[1:100],
    sequence = fasta_data$sequence[1:100],
    abundance = abundance_data$abundance[1:100],
    options = options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 10)
  expect_equal(nrow(results$chimera_report), 100)

  sequence_names <- readRDS(rchime_example("miseq_names_by_sample.rds"))
  sequences <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
  abundances <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))

  expect_error(rchimeDenovo(
    sequence_name = sequence_names[1:2],
    sequence = sequences[1:5],
    abundance = abundances
  ))

  expect_error(rchimeDenovo(
    sequence_name = sequence_names,
    sequence = sequences,
    abundance = abundances[1:2]
  ))

  options <- list(xn = 8, dereplicate = FALSE)
  results <- rchimeDenovo(
    sequence_name = sequence_names[1:2],
    sequence = sequences[1:2],
    abundance = abundances[1:2], options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 520)
  expect_equal(nrow(results$chimera_report), 966)
})

test_that("test rchimeReference ", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  results <- rchimeReference(
    sequence_name = fasta_data$sequence_name[1:100],
    sequence = fasta_data$sequence[1:100],
    reference_name = reference_data$sequence_name,
    reference_sequence = reference_data$sequence
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 24)
  expect_equal(nrow(results$chimera_report), 100)
})

test_that("test rchimeDenovoSingleSample ", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  abundance_data <- readRDS(rchime_example("single_sample_abundance.rds"))

  results <- rchimeDenovoSingleSample(
    sequence_name = fasta_data$sequence_name[100:200],
    sequence = fasta_data$sequence[100:200],
    abundance = abundance_data$abundance[100:200]
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 18)
  expect_equal(nrow(results$chimera_report), 101)
})

test_that("test rchimeDenovo ", {
  sequence_names <- readRDS(rchime_example("miseq_names_by_sample.rds"))
  sequences <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
  abundances <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))

  options <- rchime_options(dereplicate = TRUE)

  results <- rchimeDenovo(
    sequence_name = sequence_names[1:2],
    sequence = sequences[1:2],
    abundance = abundances[1:2], options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 3)
  expect_equal(length(results$chimeras), 513)
  expect_equal(nrow(results$chimera_report), 966)

  results <- rchimeDenovo(
    sequence_name = sequence_names[3:4],
    sequence = sequences[3:4],
    abundance = abundances[3:4], list(dereplicate = FALSE)
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 364)
  expect_equal(nrow(results$chimera_report), 722)
})
