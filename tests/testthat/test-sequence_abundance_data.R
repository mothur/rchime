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

  expect_equal(abunds$get_total(), 3)
  expect_equal(abunds$get_group_totals(), c(2,1))

  # set the abundance of 'seq1' in 'sample1' to 150,
  #                      'seq1' in 'sample2' to 0
  abunds$set_abundance("seq1", list(c(150, 0)))

  # set the abundance of 'seq3' in 'sample1' to 75,
  #                      'seq3' in 'sample2' to 5
  abunds$set_abundance("seq3", list(c(75, 5)))

  expect_equal(abunds$get_total(), 231)
  expect_equal(abunds$get_group_totals(), c(226,5))

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

  abunds$merge_seqs(c("seq1", "seq2"))
  expect_equal(abunds$get_total(), 231)
  expect_equal(abunds$get_group_totals(), c(226,5))
  # seq1 now represent seq1 and seq2
  expect_equal(abunds$get_abund("seq1"), 151)

  # seq1 reps seq1 and seq2
  abunds$remove_seq("seq1")
  expect_equal(abunds$get_groups(), c("sample1", "sample2"))
  expect_equal(abunds$get_num_groups(), 2)
  expect_equal(abunds$get_total(), 80)
  expect_equal(abunds$get_group_totals(), c(75,5))

  abunds$remove_seq("seq3")
  expect_equal(abunds$get_num_groups(), 0)
  expect_equal(abunds$get_total(), 0)
})
