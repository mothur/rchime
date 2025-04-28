# test sequence_abdundance_data R6 Class

test_that("test sequence_abdundance_data R6 class", {
  abunds <- sequence_abundance_data$new()
  names <- c("seq1", "seq2", "seq3")
  abunds$add_seqs(names)

  expect_equal(abunds$get_total(), 3)

  # set total abundance. seq1 -> 100, seq2 -> 50, seq3 -> 5
  new_abunds <- list(100, 50, 5)
  abunds$set_abundance("seq1", new_abunds[1])
  abunds$set_abundance("seq2", new_abunds[2])
  abunds$set_abundance("seq3", new_abunds[3])

  expect_equal(abunds$get_total(), 155)
  expect_false(abunds$has_groups())

  abunds <- sequence_abundance_data$new()
  groups <- c("sample1", "sample1", "sample2")
  abunds$set_group_assignments(names, groups)

  expect_error(abunds$add_seqs(names))
  expect_error(abunds$get_abund("bad_name"))

  abunds <- sequence_abundance_data$new()
  groups <- c("sample1", "sample1", "sample2")
  abunds$set_group_assignments(names, groups)

  expect_equal(abunds$get_total(), 3)
  expect_equal(abunds$get_group_totals(), c(2, 1))

  # set the abundance of 'seq1' in 'sample1' to 150,
  #                      'seq1' in 'sample2' to 0
  abunds$set_abundance("seq1", list(c(150, 0)))

  # set the abundance of 'seq3' in 'sample1' to 75,
  #                      'seq3' in 'sample2' to 5
  abunds$set_abundance("seq3", list(c(75, 5)))

  expect_equal(abunds$get_total(), 231)
  expect_equal(abunds$get_group_totals(), c(226, 5))

  expect_equal(abunds$get_abund("seq1"), 150)
  expect_equal(abunds$get_abund("seq2"), 1)
  expect_equal(abunds$get_abund("seq3"), 80)

  expect_equal(abunds$get_abunds("seq1"), c(150, 0))
  expect_equal(abunds$get_abunds("seq2"), c(1, 0))
  expect_equal(abunds$get_abunds("seq3"), c(75, 5))

  expect_equal(abunds$get_groups(), c("sample1", "sample2"))
  expect_equal(abunds$get_groups("seq1"), c("sample1"))
  expect_equal(abunds$get_groups("seq2"), c("sample1"))
  expect_equal(abunds$get_groups("seq3"), c("sample1", "sample2"))
  expect_equal(abunds$get_num_groups(), 2)

  abunds_by_sample <- abunds$get_seqs_abunds(names, bysample = TRUE)
  # seq1 abunds by sample
  expect_equal(abunds_by_sample[[1]], c(150, 0))
  # seq2 abunds by sample
  expect_equal(abunds_by_sample[[2]], c(1, 0))
  # seq3 abunds by sample
  expect_equal(abunds_by_sample[[3]], c(75, 5))
  expect_equal(abunds$get_seqs_abunds(names), c(150, 1, 80))

  expect_true(abunds$has_groups())
  expect_false(abunds$has_group("bad_group"))
  expect_false(abunds$has_group("sample2", "seq2"))

  # merge sample1
  abunds$merge_seqs(c("seq1", "seq2"), group = "sample1")
  expect_equal(abunds$get_total(), 231)
  expect_equal(abunds$get_group_totals(), c(226, 5))
  # seq1 now represent seq1 and seq2
  expect_equal(abunds$get_abunds("seq1"), c(151, 0))

  abunds$merge_seqs(c("seq1", "seq2"))
  expect_equal(abunds$get_total(), 231)
  expect_equal(abunds$get_group_totals(), c(226, 5))
  # seq1 now represent seq1 and seq2
  expect_equal(abunds$get_abund("seq1"), 151)

  # seq1 reps seq1 and seq2
  abunds$remove_seq("seq1")
  expect_equal(abunds$get_groups(), c("sample1", "sample2"))
  expect_equal(abunds$get_num_groups(), 2)
  expect_equal(abunds$get_total(), 80)
  expect_equal(abunds$get_group_totals(), c(75, 5))

  abunds$remove_seq("seq3")
  expect_equal(abunds$get_num_groups(), 0)
  expect_equal(abunds$get_total(), 0)
})

test_that("test set_group_assignments - names, groups, abundances", {
  # mothur count file
  # "Representative_Sequence     total   sample2	sample3	sample4",
  # "seq1	1150	250	400	500", "seq2	115	25	40	50",
  # "seq3	50	25	25	0", "seq4	4	0	0	4"

  # inputs as sample table
  names <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3",
    "seq4"
  )
  groups <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3",
    "sample4"
  )
  abundances <- c(
    250, 400, 500,
    25, 40, 50,
    25, 25,
    4
  )

  abunds <- sequence_abundance_data$new()
  abunds$set_group_assignments(names, groups, abundances)

  expect_equal(abunds$get_total(), 1319)
  expect_equal(abunds$get_num_groups(), 3)

  group_totals <- abunds$get_group_totals()

  expect_equal(group_totals[1], 300)
  expect_equal(group_totals[2], 465)
  expect_equal(group_totals[3], 554)

  sample_table <- abunds$get_sample_table()

  expect_equal(sample_table$name, names)
  expect_equal(sample_table$group, groups)
  expect_equal(sample_table$abundance, abundances)

  abunds <- sequence_abundance_data$new()
  abunds$set_group_assignments(names, groups)

  expect_equal(abunds$get_total(), 9)
  expect_equal(abunds$get_num_groups(), 3)

  group_totals <- abunds$get_group_totals()

  expect_equal(group_totals[1], 3)
  expect_equal(group_totals[2], 3)
  expect_equal(group_totals[3], 3)

  abunds <- sequence_abundance_data$new()
  abunds$add_seqs(unique(names))
  abunds$set_abundance("seq1", list(c(100)))
  abunds$set_abundance("seq2", list(c(200)))
  abunds$set_abundance("seq3", list(c(300)))
  abunds$set_abundance("seq4", list(c(400)))

  expect_equal(abunds$get_total(), 1000)
  expect_equal(abunds$get_num_groups(), 0)

  sample_table <- abunds$get_sample_table()

  expect_equal(sample_table$name, unique(names))
  expect_equal(sample_table$abundance, c(100, 200, 300, 400))
})

test_that("test set_group_assignments - export", {
  # no groups
  abunds <- sequence_abundance_data$new()
  abunds$set_group_assignments(
    filename =
      rchime_example("test_nogroups.count_table")
  )

  df <- abunds$export()

  expect_equal(colnames(df), c("total_abundance"))
  expect_equal(df[[1]], c(108939, 1734, 4773, 15750, 150122, 14378))

  # groups
  # inputs as sample table
  names <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3",
    "seq4"
  )
  groups <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3",
    "sample4"
  )
  abundances <- c(
    250, 400, 500,
    25, 40, 50,
    25, 25,
    4
  )

  abunds <- sequence_abundance_data$new()
  abunds$set_group_assignments(names, groups, abundances)

  df <- abunds$export()

  expect_equal(abunds$get_total(), 1319)
  expect_equal(abunds$get_num_groups(), 3)

  expect_equal(colnames(df), c(
    "total_abundance", "sample2", "sample3",
    "sample4"
  ))

  expect_equal(df[[1]], c(1150, 115, 50, 4))
  expect_equal(df[[2]], c(250, 25, 25, 0))
  expect_equal(df[[3]], c(400, 40, 25, 0))
  expect_equal(df[[4]], c(500, 50, 0, 4))
})
