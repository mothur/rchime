# test rpp_rchime

test_that("test rchimeReference - errors", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  # not equal sizes
  expect_error(rchimeReference(
    sequence_names = fasta_data$sequence_names[1:50],
    sequences = fasta_data$sequences[1:100],
    reference_names = reference_data$sequence_names,
    reference_sequences = reference_data$sequences
  ))

  # not equal sizes
  expect_error(rchimeReference(
    sequence_names = fasta_data$sequence_names[1:100],
    sequences = fasta_data$sequences[1:100],
    reference_names = reference_data$sequence_names[1:5],
    reference_sequences = reference_data$sequences[1:10]
  ))

  options <- list(xn = 8.0)
  # create options but leave off processors and dereplicate
  results <- rchimeReference(
    sequence_names = fasta_data$sequence_names[1:100],
    sequences = fasta_data$sequences[1:100],
    reference_names = reference_data$sequence_names,
    reference_sequences = reference_data$sequences,
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
    sequence_names = fasta_data$sequence_names[1:50],
    sequences = fasta_data$sequences[1:100],
    abundances = abundance_data$abundances
  ))

  # not equal sizes
  expect_error(rchimeDenovoSingleSample(
    sequence_names = fasta_data$sequence_names[1:100],
    sequences = fasta_data$sequences[1:100],
    abundances = abundance_data$abundances
  ))

  options <- list(xn = 8.0)
  # create options but leave off processors and dereplicate
  results <- rchimeDenovoSingleSample(
    sequence_names = fasta_data$sequence_names,
    sequences = fasta_data$sequences,
    abundances = abundance_data$abundances,
    options = options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 3719)
  expect_equal(nrow(results$chimera_report), 6084)

  sequence_names <- readRDS(rchime_example("miseq_names_by_sample.rds"))
  sequences <- readRDS(rchime_example("miseq_sequences_by_sample.rds"))
  abundances <- readRDS(rchime_example("miseq_abundance_by_sample.rds"))

  expect_error(rchimeDenovo(
    sequence_names = sequence_names[1:2],
    sequences = sequences[1:5],
    abundances = abundances
  ))

  expect_error(rchimeDenovo(
    sequence_names = sequence_names,
    sequences = sequences,
    abundances = abundances[1:2]
  ))

  options <- list(xn = 8)
  results <- rchimeDenovo(
    sequence_names = sequence_names,
    sequences = sequences,
    abundances = abundances, options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 3751)
  expect_equal(nrow(results$chimera_report), 6084)
})

test_that("test rchimeReference ", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  results <- rchimeReference(
    sequence_names = fasta_data$sequence_names[1:100],
    sequences = fasta_data$sequences[1:100],
    reference_names = reference_data$sequence_names,
    reference_sequences = reference_data$sequences
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
    sequence_names = fasta_data$sequence_names[100:200],
    sequences = fasta_data$sequences[100:200],
    abundances = abundance_data$abundances[100:200]
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
    sequence_names = sequence_names[1:2],
    sequences = sequences[1:2],
    abundances = abundances[1:2], options
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 3)
  expect_equal(length(results$chimeras), 513)
  expect_equal(nrow(results$chimera_report), 966)

  results <- rchimeDenovo(
    sequence_names = sequence_names[3:10],
    sequences = sequences[3:10],
    abundances = abundances[3:10]
  )

  # checks to make sure the correct things are created
  expect_equal(length(results), 2)
  expect_equal(length(results$chimeras), 1648)
  expect_equal(nrow(results$chimera_report), 2901)
})
