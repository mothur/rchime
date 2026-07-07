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
  reference_data <- silva_gold()

  chimera_report <- rchime(
    data = fasta_data[1:100, ],
    reference = reference_data,
    rchime_options = options,
    verbose = FALSE
  )

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 2)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime by reference - data.frame", {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  reference_data <- silva_gold()

  chimera_report <- rchime(
    data = fasta_data[1:100, ],
    reference = reference_data
  )

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 12)
  expect_equal(nrow(chimera_report$chimera_report), 100)

  # spot check chimera report
  # check first chimeric sequence
  first_chimeras_name <- "M00967_43_000000000-A3JHG_1_2107_6538_14332"

  # query and parent names
  expect_equal(chimera_report$chimeras[1], first_chimeras_name)
  expect_equal(chimera_report$chimera_report[[14, 2]], first_chimeras_name)
  expect_equal(chimera_report$chimera_report[[14, 3]], "S000013923")
  expect_equal(chimera_report$chimera_report[[14, 4]], "S000022285")
  expect_equal(chimera_report$chimera_report[[14, 5]], "S000022285")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[14, 11]], 22)
  expect_equal(chimera_report$chimera_report[[14, 12]], 6)
  expect_equal(chimera_report$chimera_report[[14, 13]], 14)
  expect_equal(chimera_report$chimera_report[[14, 14]], 19)
  expect_equal(chimera_report$chimera_report[[14, 15]], 0)
  expect_equal(chimera_report$chimera_report[[14, 16]], 5)

  # check first non chimeric sequence
  first_non_chimeras_name <- "M00967_43_000000000-A3JHG_1_1103_18369_27790"

  expect_equal(chimera_report$chimera_report[[2, 2]], first_non_chimeras_name)
  expect_equal(chimera_report$chimera_report[[2, 3]], "*")
  expect_equal(chimera_report$chimera_report[[2, 4]], "*")
  expect_equal(chimera_report$chimera_report[[2, 5]], "*")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[2, 11]], -1)
  expect_equal(chimera_report$chimera_report[[2, 12]], -1)
  expect_equal(chimera_report$chimera_report[[2, 13]], -1)
  expect_equal(chimera_report$chimera_report[[2, 14]], -1)
  expect_equal(chimera_report$chimera_report[[2, 15]], -1)
  expect_equal(chimera_report$chimera_report[[2, 16]], -1)

  # first non chimera with parents selected
  first_non_chimera_wp <- "M00967_43_000000000-A3JHG_1_1109_8934_16443"

  expect_equal(chimera_report$chimera_report[[3, 2]], first_non_chimera_wp)
  expect_equal(chimera_report$chimera_report[[3, 3]], "S000013923")
  expect_equal(chimera_report$chimera_report[[3, 4]], "7000004131498097")
  expect_equal(chimera_report$chimera_report[[3, 5]], "7000004131498097")
  expect_equal(chimera_report$chimera_report[[3, 11]], 5)
  expect_equal(chimera_report$chimera_report[[3, 12]], 0)
  expect_equal(chimera_report$chimera_report[[3, 13]], 4)
  expect_equal(chimera_report$chimera_report[[3, 14]], 13)
  expect_equal(chimera_report$chimera_report[[3, 15]], 3)
  expect_equal(chimera_report$chimera_report[[3, 16]], 23)
})

test_that("test rchime de novo - data.frame - single sample ", {
  # no need to check dereplicate since one sample

  data <- readRDS(rchime_example("miseq_data_frame.rds"))

  chimera_report <- rchime(data[1:100, ], verbose = FALSE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 19)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime de novo -data.frame/mult samples, dereplicate = TRUE", {
  data <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))
  chimera_report <- rchime(data, dereplicate = TRUE, verbose = FALSE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 3)
  expect_equal(length(chimera_report$chimeras), 90)
  expect_equal(nrow(chimera_report$chimera_report), 500)
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
  chimeras_name <- "M00967_43_000000000-A3JHG_1_1101_10331_23332"

  # query and parent names
  expect_equal(chimera_report$chimera_report[[20, 2]], chimeras_name)
  expect_equal(
    chimera_report$chimera_report[[20, 3]],
    "M00967_43_000000000-A3JHG_1_1103_13364_5496"
  )
  expect_equal(
    chimera_report$chimera_report[[20, 4]],
    "M00967_43_000000000-A3JHG_1_1101_23238_24359"
  )
  expect_equal(
    chimera_report$chimera_report[[20, 5]],
    "M00967_43_000000000-A3JHG_1_1101_23238_24359"
  )

  # right and left votes
  expect_equal(chimera_report$chimera_report[[20, 11]], 8)
  expect_equal(chimera_report$chimera_report[[20, 12]], 0)
  expect_equal(chimera_report$chimera_report[[20, 13]], 0)
  expect_equal(chimera_report$chimera_report[[20, 14]], 49)
  expect_equal(chimera_report$chimera_report[[20, 15]], 6)
  expect_equal(chimera_report$chimera_report[[20, 16]], 7)

  # check first non chimeric sequence
  first_non_chimeras_name <- "M00967_43_000000000-A3JHG_1_1101_18089_2781"

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
  first_non_chimera_wp <- "M00967_43_000000000-A3JHG_1_1103_19870_21567"

  expect_equal(chimera_report$chimera_report[[7, 2]], first_non_chimera_wp)
  expect_equal(
    chimera_report$chimera_report[[7, 3]],
    "M00967_43_000000000-A3JHG_1_1101_18922_4934"
  )
  expect_equal(
    chimera_report$chimera_report[[7, 4]],
    "M00967_43_000000000-A3JHG_1_1101_18089_2781"
  )
  expect_equal(
    chimera_report$chimera_report[[7, 5]],
    "M00967_43_000000000-A3JHG_1_1101_18922_4934"
  )
  expect_equal(chimera_report$chimera_report[[7, 11]], 6)
  expect_equal(chimera_report$chimera_report[[7, 12]], 2)
  expect_equal(chimera_report$chimera_report[[7, 13]], 14)
  expect_equal(chimera_report$chimera_report[[7, 14]], 3)
  expect_equal(chimera_report$chimera_report[[7, 15]], 2)
  expect_equal(chimera_report$chimera_report[[7, 16]], 5)
})

test_that("test rchime de novo strollur/mulitple samples, dereplicate = FALSE", {
  data <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))

  options <- rchime_options(dereplicate = FALSE)
  chimera_report <- rchime(data, verbose = FALSE, rchime_options = options)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 108)
  expect_equal(nrow(chimera_report$chimera_report), 500)

  # spot check chimera report
  # check first chimeric sequence
  chimeras_name <- "M00967_43_000000000-A3JHG_1_1101_10331_23332"

  # query and parent names
  expect_equal(chimera_report$chimeras[1], chimeras_name)
  expect_equal(chimera_report$chimera_report[[18, 2]], chimeras_name)
  expect_equal(
    chimera_report$chimera_report[[18, 3]],
    "M00967_43_000000000-A3JHG_1_1103_13364_5496"
  )
  expect_equal(
    chimera_report$chimera_report[[18, 4]],
    "M00967_43_000000000-A3JHG_1_1101_23238_24359"
  )
  expect_equal(
    chimera_report$chimera_report[[18, 5]],
    "M00967_43_000000000-A3JHG_1_1101_23238_24359"
  )

  # right and left votes
  expect_equal(chimera_report$chimera_report[[18, 11]], 8)
  expect_equal(chimera_report$chimera_report[[18, 12]], 0)
  expect_equal(chimera_report$chimera_report[[18, 13]], 0)
  expect_equal(chimera_report$chimera_report[[18, 14]], 49)
  expect_equal(chimera_report$chimera_report[[18, 15]], 6)
  expect_equal(chimera_report$chimera_report[[18, 16]], 7)

  # check first non chimeric sequence
  first_non_chimeras_name <- "M00967_43_000000000-A3JHG_1_1101_18089_2781"

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

  # first non chimera with parents selected
  first_non_chimera_wp <- "M00967_43_000000000-A3JHG_1_1103_19870_21567"

  expect_equal(chimera_report$chimera_report[[7, 2]], first_non_chimera_wp)
  expect_equal(
    chimera_report$chimera_report[[7, 3]],
    "M00967_43_000000000-A3JHG_1_1101_18922_4934"
  )
  expect_equal(
    chimera_report$chimera_report[[7, 4]],
    "M00967_43_000000000-A3JHG_1_1101_18089_2781"
  )
  expect_equal(
    chimera_report$chimera_report[[7, 5]],
    "M00967_43_000000000-A3JHG_1_1101_18922_4934"
  )
  expect_equal(chimera_report$chimera_report[[7, 11]], 6)
  expect_equal(chimera_report$chimera_report[[7, 12]], 2)
  expect_equal(chimera_report$chimera_report[[7, 13]], 14)
  expect_equal(chimera_report$chimera_report[[7, 14]], 3)
  expect_equal(chimera_report$chimera_report[[7, 15]], 2)
  expect_equal(chimera_report$chimera_report[[7, 16]], 5)
})
