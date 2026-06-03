# test rchime function

test_that("test rchime - errors and rchime_options", {
  data <- strollur::load_dataset(
    rchime_example("strollur_miseq_tiny.rds")
  )
  reference_ob <- strollur::load_dataset(
    rchime_example("strollur_reference.rds")
  )

  fasta_data <- strollur::report(data, type = "fasta")

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
  strollur::add(data, table = fasta_data[1:3, ], type = "sequence")

  chimera_report <- rchime(data, reference = reference_ob)

  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 0)
  expect_equal(nrow(chimera_report$chimera_report), 3)

  options <- rchime_options(
    minh = 0.6, abskew = 1.9, mindiv = 0.505,
    xn = 8.01, dn = 1.401, maxp = 4
  )

  data <- strollur::load_dataset(rchime_example("strollur_miseq_tiny.rds"))

  chimera_report <- rchime(data,
    reference = reference_ob,
    rchime_options = options,
    silent = TRUE
  )

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 13)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime by reference - strollur", {
  data <- strollur::load_dataset(rchime_example("strollur_miseq_tiny.rds"))
  reference <- strollur::load_dataset(rchime_example("strollur_reference.rds"))

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
  # check last chimeric sequence
  last_chimeras_name <- "M00967_43_000000000-A3JHG_1_2113_29036_16812"

  # query and parent names
  expect_equal(chimera_report$chimeras[24], last_chimeras_name)
  expect_equal(chimera_report$chimera_report[[99, 2]], last_chimeras_name)
  expect_equal(chimera_report$chimera_report[[99, 3]], "AJ404681.1")
  expect_equal(chimera_report$chimera_report[[99, 4]], "AJ400237.1")
  expect_equal(chimera_report$chimera_report[[99, 5]], "AJ404681.1")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[99, 11]], 22)
  expect_equal(chimera_report$chimera_report[[99, 12]], 2)
  expect_equal(chimera_report$chimera_report[[99, 13]], 11)
  expect_equal(chimera_report$chimera_report[[99, 14]], 7)
  expect_equal(chimera_report$chimera_report[[99, 15]], 0)
  expect_equal(chimera_report$chimera_report[[99, 16]], 0)

  # check non chimeric sequence
  non_chimeras_name <- "M00967_43_000000000-A3JHG_1_2106_10410_19621"

  expect_equal(chimera_report$chimera_report[[75, 2]], non_chimeras_name)
  expect_equal(chimera_report$chimera_report[[75, 3]], "*")
  expect_equal(chimera_report$chimera_report[[75, 4]], "*")
  expect_equal(chimera_report$chimera_report[[75, 5]], "*")

  # right and left votes
  expect_equal(chimera_report$chimera_report[[75, 11]], -1)
  expect_equal(chimera_report$chimera_report[[75, 12]], -1)
  expect_equal(chimera_report$chimera_report[[75, 13]], -1)
  expect_equal(chimera_report$chimera_report[[75, 14]], -1)
  expect_equal(chimera_report$chimera_report[[75, 15]], -1)
  expect_equal(chimera_report$chimera_report[[75, 16]], -1)

  # non chimera with parents selected
  non_chimera_wp <- "M00967_43_000000000-A3JHG_1_1109_8934_16443"

  expect_equal(chimera_report$chimera_report[[28, 2]], non_chimera_wp)
  expect_equal(chimera_report$chimera_report[[28, 3]], "AJ400264.1")
  expect_equal(chimera_report$chimera_report[[28, 4]], "AJ400267.1")
  expect_equal(chimera_report$chimera_report[[28, 5]], "AJ400264.1")
  expect_equal(chimera_report$chimera_report[[28, 11]], 10)
  expect_equal(chimera_report$chimera_report[[28, 12]], 0)
  expect_equal(chimera_report$chimera_report[[28, 13]], 1)
  expect_equal(chimera_report$chimera_report[[28, 14]], 7)
  expect_equal(chimera_report$chimera_report[[28, 15]], 1)
  expect_equal(chimera_report$chimera_report[[28, 16]], 5)

  # check by reference with multiple sample object
  data <- strollur::load_dataset(
    rchime_example("strollur_multi_sample_tiny.rds")
  )
  chimera_report <- rchime(data, reference = reference)

  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 24)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime denovo - strollur - single sample ", {
  data <- strollur::load_dataset(
    rchime_example("strollur_single_sample_tiny.rds")
  )

  chimera_report <- rchime(data, silent = TRUE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 10)
  expect_equal(nrow(chimera_report$chimera_report), 100)
})

test_that("test rchime denovo -strollur/mulitple samples, dereplicate = TRUE", {
  data <- strollur::load_dataset(
    rchime_example("strollur_multi_sample_small.rds")
  )
  chimera_report <- rchime(data, silent = TRUE)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
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

  # spot check chimera report
  # check chimeric sequence
  chimeras_name <- "M00967_43_000000000-A3JHG_1_1101_10331_23332"

  # query and parent names
  expect_equal(chimera_report$chimeras[1], chimeras_name)
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

test_that("test rchime denovo strollur/mulitple samples, dereplicate = FALSE", {
  data <- strollur::load_dataset(
    rchime_example("strollur_multi_sample_small.rds")
  )

  options <- rchime_options(dereplicate = FALSE)
  chimera_report <- rchime(data, silent = TRUE, rchime_options = options)

  # checks to make sure the correct things are created
  expect_equal(length(chimera_report), 2)
  expect_equal(length(chimera_report$chimeras), 108)
  expect_equal(nrow(chimera_report$chimera_report), 500)

  # spot check chimera report
  # check chimeric sequence
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
