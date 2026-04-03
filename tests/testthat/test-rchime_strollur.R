# test rchime function

test_that("test rchime - errors and rchime_options", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  data <- strollur::new_dataset("rchime reference example")
  strollur::add(data, table = fasta_data[1:100, ], type = "sequences")

  reference_ob <- new_dataset("Silva V4 Region")
  strollur::add(reference_ob, table = reference_data, type = "sequences")

  # this should not change since the rchime command below does not remove
  num_seqs <- strollur::count(data)

  options <- rchime_options(dereplicate = FALSE)
  chimera_report <- rchime(data,
    reference = reference_ob,
    remove_chimeras = FALSE,
    silent = FALSE,
    rchime_options = options
  )

  expect_equal(strollur::count(data), num_seqs)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 24)
  expect_equal(nrow(chimera_report$chimera_report), 100)

  expect_error(rchime(data, reference = "not a strollur object"))

  # test small enough set no chimeras are detected
  data <- strollur::new_dataset("rchime reference example")
  strollur::add(data, table = fasta_data[1:3, ], type = "sequences")

  chimera_report <- rchime(data, reference = reference_ob)

  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 0)
  expect_equal(nrow(chimera_report$chimera_report), 3)

  options <- rchime_options(
    minh = 0.6, abskew = 1.9, mindiv = 0.505,
    xn = 8.01, dn = 1.401, maxp = 4
  )

  data <- strollur::new_dataset("rchime reference example")
  strollur::add(data, table = fasta_data[1:100, ], type = "sequences")

  chimera_report <- rchime(data,
    reference = reference_ob,
    rchime_options = options,
    silent = TRUE
  )

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 14)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime by reference - strollur", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  data <- strollur::new_dataset("rchime reference example")
  strollur::add(data, table = fasta_data[1:100, ], type = "sequences")

  reference <- new_dataset("Silva V4 Region")
  strollur::add(reference, table = reference_data, type = "sequences")

  chimera_report <- rchime(data, reference = reference, remove_chimeras = FALSE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 24)
  expect_equal(nrow(chimera_report$chimera_report), 100)

  chimera_report <- rchime(data, reference = reference)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 24)
  expect_equal(nrow(chimera_report$chimera_report), 100)


  # spot check chimera report
  # check first chimeric sequence
  first_chimeras_name <- "M00967_43_000000000-A3JHG_1_2113_29036_16812"

  # query and parent names
  expect_equal(chimera_report$chimeras[1], first_chimeras_name)
  expect_equal(chimera_report$chimera_report[[4, 2]], first_chimeras_name)
  expect_equal(chimera_report$chimera_report[[4, 3]], "AJ404681.1")
  expect_equal(chimera_report$chimera_report[[4, 4]], "AJ400237.1")
  expect_equal(chimera_report$chimera_report[[4, 5]], "AJ404681.1")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[4, 11]], 22)
  expect_equal(chimera_report$chimera_report[[4, 12]], 2)
  expect_equal(chimera_report$chimera_report[[4, 13]], 11)
  expect_equal(chimera_report$chimera_report[[4, 14]], 7)
  expect_equal(chimera_report$chimera_report[[4, 15]], 0)
  expect_equal(chimera_report$chimera_report[[4, 16]], 0)

  # check first non chimeric sequence
  first_non_chimeras_name <- "M00967_43_000000000-A3JHG_1_2106_10410_19621"

  expect_equal(chimera_report$chimera_report[[1, 2]], first_non_chimeras_name)
  expect_equal(chimera_report$chimera_report[[1, 3]], "*")
  expect_equal(chimera_report$chimera_report[[1, 4]], "*")
  expect_equal(chimera_report$chimera_report[[1, 5]], "*")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[1, 11]], -1)
  expect_equal(chimera_report$chimera_report[[1, 12]], -1)
  expect_equal(chimera_report$chimera_report[[1, 13]], -1)
  expect_equal(chimera_report$chimera_report[[1, 14]], -1)
  expect_equal(chimera_report$chimera_report[[1, 15]], -1)
  expect_equal(chimera_report$chimera_report[[1, 16]], -1)

  # first non chimera with parents selected
  first_non_chimera_wp <- "M00967_43_000000000-A3JHG_1_1109_8934_16443"

  expect_equal(chimera_report$chimera_report[[3, 2]], first_non_chimera_wp)
  expect_equal(chimera_report$chimera_report[[3, 3]], "AJ400264.1")
  expect_equal(chimera_report$chimera_report[[3, 4]], "AJ400267.1")
  expect_equal(chimera_report$chimera_report[[3, 5]], "AJ400264.1")
  expect_equal(chimera_report$chimera_report[[3, 11]], 10)
  expect_equal(chimera_report$chimera_report[[3, 12]], 0)
  expect_equal(chimera_report$chimera_report[[3, 13]], 1)
  expect_equal(chimera_report$chimera_report[[3, 14]], 7)
  expect_equal(chimera_report$chimera_report[[3, 15]], 1)
  expect_equal(chimera_report$chimera_report[[3, 16]], 5)
})

test_that("test rchime denovo - strollur - single sample ", {
  # no need to check dereplicate since one sample

  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  abundance_data <- readRDS(rchime_example("single_sample_abundance.rds"))

  data <- strollur::new_dataset("rchime denovo example")
  strollur::add(data, table = fasta_data, type = "sequences")
  strollur::assign(data, table = abundance_data, type = "sequence_abundance")

  chimera_report <- rchime(data, silent = TRUE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 3719)
  expect_equal(nrow(chimera_report$chimera_report), 6084)
})

test_that("test rchime denovo -strollur/mulitple samples, dereplicate = TRUE", {
  data <- strollur::load_dataset(rchime_example("strollur_multi_sample.rds"))
  chimera_report <- rchime(data, silent = TRUE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 3588)
  expect_equal(nrow(chimera_report$chimera_report), 6084)
  expect_equal(
    names(chimera_report$chimera_report),
    c(
      "Score", "Query", "ParentA", "ParentB", "Top_Parent",
      "QM", "QA", "QB", "QAB", "QT",
      "LY", "LN", "LA", "RY", "RN", "RA", "Div", "Chimeric_Status"
    )
  )

  # spot check chimera report
  # check first chimeric sequence
  first_chimeras_name <- "M00967_43_000000000-A3JHG_1_1106_11629_14238"

  # query and parent names
  expect_equal(chimera_report$chimeras[1], first_chimeras_name)
  expect_equal(chimera_report$chimera_report[[66, 2]], first_chimeras_name)
  expect_equal(
    chimera_report$chimera_report[[66, 3]],
    "M00967_43_000000000-A3JHG_1_1107_15750_18592"
  )
  expect_equal(
    chimera_report$chimera_report[[66, 4]],
    "M00967_43_000000000-A3JHG_1_2101_22400_13416"
  )
  expect_equal(
    chimera_report$chimera_report[[66, 5]],
    "M00967_43_000000000-A3JHG_1_1107_15750_18592"
  )

  # right and left votes
  expect_equal(chimera_report$chimera_report[[66, 11]], 13)
  expect_equal(chimera_report$chimera_report[[66, 12]], 0)
  expect_equal(chimera_report$chimera_report[[66, 13]], 0)
  expect_equal(chimera_report$chimera_report[[66, 14]], 4)
  expect_equal(chimera_report$chimera_report[[66, 15]], 0)
  expect_equal(chimera_report$chimera_report[[66, 16]], 1)

  # check first non chimeric sequence
  first_non_chimeras_name <- "M00967_43_000000000-A3JHG_1_1107_15750_18592"

  expect_equal(chimera_report$chimera_report[[1, 2]], first_non_chimeras_name)
  expect_equal(chimera_report$chimera_report[[1, 3]], "*")
  expect_equal(chimera_report$chimera_report[[1, 4]], "*")
  expect_equal(chimera_report$chimera_report[[1, 5]], "*")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[1, 11]], -1)
  expect_equal(chimera_report$chimera_report[[1, 12]], -1)
  expect_equal(chimera_report$chimera_report[[1, 13]], -1)
  expect_equal(chimera_report$chimera_report[[1, 14]], -1)
  expect_equal(chimera_report$chimera_report[[1, 15]], -1)
  expect_equal(chimera_report$chimera_report[[1, 16]], -1)

  # first non chimera with parents selected
  first_non_chimera_wp <- "M00967_43_000000000-A3JHG_1_2108_21552_11279"

  expect_equal(chimera_report$chimera_report[[8, 2]], first_non_chimera_wp)
  expect_equal(
    chimera_report$chimera_report[[8, 3]],
    "M00967_43_000000000-A3JHG_1_2110_12856_16229"
  )
  expect_equal(
    chimera_report$chimera_report[[8, 4]],
    "M00967_43_000000000-A3JHG_1_1112_6862_18037"
  )
  expect_equal(
    chimera_report$chimera_report[[8, 5]],
    "M00967_43_000000000-A3JHG_1_2110_12856_16229"
  )
  expect_equal(chimera_report$chimera_report[[8, 11]], 8)
  expect_equal(chimera_report$chimera_report[[8, 12]], 2)
  expect_equal(chimera_report$chimera_report[[8, 13]], 17)
  expect_equal(chimera_report$chimera_report[[8, 14]], 5)
  expect_equal(chimera_report$chimera_report[[8, 15]], 1)
  expect_equal(chimera_report$chimera_report[[8, 16]], 21)
})

test_that("test rchime denovo strollur/mulitple samples, dereplicate = FALSE", {
  data <- strollur::load_dataset(rchime_example("strollur_multi_sample.rds"))

  options <- rchime_options(dereplicate = FALSE)
  chimera_report <- rchime(data, silent = TRUE, rchime_options = options)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 3751)
  expect_equal(nrow(chimera_report$chimera_report), 6084)

  # spot check chimera report
  # check first chimeric sequence
  first_chimeras_name <- "M00967_43_000000000-A3JHG_1_2104_18459_1792"

  # query and parent names
  expect_equal(chimera_report$chimeras[1], first_chimeras_name)
  expect_equal(chimera_report$chimera_report[[55, 2]], first_chimeras_name)
  expect_equal(
    chimera_report$chimera_report[[55, 3]],
    "M00967_43_000000000-A3JHG_1_1105_17475_15394"
  )
  expect_equal(
    chimera_report$chimera_report[[55, 4]],
    "M00967_43_000000000-A3JHG_1_1113_18037_24127"
  )
  expect_equal(
    chimera_report$chimera_report[[55, 5]],
    "M00967_43_000000000-A3JHG_1_1105_17475_15394"
  )

  # right and left votes
  expect_equal(chimera_report$chimera_report[[55, 11]], 28)
  expect_equal(chimera_report$chimera_report[[55, 12]], 0)
  expect_equal(chimera_report$chimera_report[[55, 13]], 0)
  expect_equal(chimera_report$chimera_report[[55, 14]], 11)
  expect_equal(chimera_report$chimera_report[[55, 15]], 3)
  expect_equal(chimera_report$chimera_report[[55, 16]], 2)

  # check first non chimeric sequence
  first_non_chimeras_name <- "M00967_43_000000000-A3JHG_1_1107_15750_18592"

  expect_equal(chimera_report$chimera_report[[1, 2]], first_non_chimeras_name)
  expect_equal(chimera_report$chimera_report[[1, 3]], "*")
  expect_equal(chimera_report$chimera_report[[1, 4]], "*")
  expect_equal(chimera_report$chimera_report[[1, 5]], "*")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[1, 11]], -1)
  expect_equal(chimera_report$chimera_report[[1, 12]], -1)
  expect_equal(chimera_report$chimera_report[[1, 13]], -1)
  expect_equal(chimera_report$chimera_report[[1, 14]], -1)
  expect_equal(chimera_report$chimera_report[[1, 15]], -1)
  expect_equal(chimera_report$chimera_report[[1, 16]], -1)

  # first non chimera with parents selected
  first_non_chimera_wp <- "M00967_43_000000000-A3JHG_1_2108_21552_11279"

  expect_equal(chimera_report$chimera_report[[8, 2]], first_non_chimera_wp)
  expect_equal(
    chimera_report$chimera_report[[8, 3]],
    "M00967_43_000000000-A3JHG_1_2110_12856_16229"
  )
  expect_equal(
    chimera_report$chimera_report[[8, 4]],
    "M00967_43_000000000-A3JHG_1_1112_6862_18037"
  )
  expect_equal(
    chimera_report$chimera_report[[8, 5]],
    "M00967_43_000000000-A3JHG_1_2110_12856_16229"
  )
  expect_equal(chimera_report$chimera_report[[8, 11]], 8)
  expect_equal(chimera_report$chimera_report[[8, 12]], 2)
  expect_equal(chimera_report$chimera_report[[8, 13]], 17)
  expect_equal(chimera_report$chimera_report[[8, 14]], 5)
  expect_equal(chimera_report$chimera_report[[8, 15]], 1)
  expect_equal(chimera_report$chimera_report[[8, 16]], 21)
})
