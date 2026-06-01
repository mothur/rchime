# test rchime function

test_that("test rchime - errors", {
  expect_error(rchime(data, reference = "not a data.frame"))
  expect_error(rchime(
    data = "not a data.frame",
    reference = "not a data.frame"
  ))

  df <- data.frame(
    not_right_tag = c("seq1", "seq2"),
    sequence = c("ATGC", "ATGC"),
    abundance = c(100, 10)
  )

  expect_error(rchime(data = df))

  options <- rchime_options(
    minh = 0.6, abskew = 1.9, mindiv = 0.505,
    xn = 8.01, dn = 1.401, maxp = 4
  )

  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  chimera_report <- rchime(
    data = fasta_data[1:100, ],
    reference = reference_data,
    rchime_options = options,
    silent = TRUE
  )

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 13)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime by reference - data.frame", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- readRDS(rchime_example("reference.rds"))

  chimera_report <- rchime(
    data = fasta_data[1:100, ],
    reference = reference_data
  )

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

  data <- readRDS(rchime_example("miseq_data_frame.rds"))

  chimera_report <- rchime(data, silent = TRUE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 3719)
  expect_equal(nrow(chimera_report$chimera_report), 6084)
})

test_that("test rchime denovo -strollur/mulitple samples, dereplicate = TRUE", {
  data <- readRDS(rchime_example("miseq_data_frame_by_sample.rds"))
  chimera_report <- rchime(data, dereplicate = TRUE, silent = TRUE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 3)
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

  samples <- c(
    "F3D0", "F3D1", "F3D141", "F3D142", "F3D143", "F3D144",
    "F3D145", "F3D146", "F3D147", "F3D148", "F3D149", "F3D150",
    "F3D2", "F3D3", "F3D5", "F3D6", "F3D7", "F3D8", "F3D9", "Mock"
  )
  expect_equal(chimera_report$set_abundance_values$samples, samples)

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
  data <- readRDS(rchime_example("miseq_data_frame_by_sample.rds"))

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
