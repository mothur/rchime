# test sequence_data_vector R6 Class

test_that("test sequence_data_vector R6 class", {
  expect_error(dataset <- sequence_data_vector$new(sequences = "seq1"))

  sequences <- c(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    ".......GAGC-TTATCCGGATT-ATTG-GTT-AAA----"
  )
  names <- c("mySeq", "seq2", "seq3")
  comments <- c("This is my comment", "", "")

  bad_data <- data.table(
    bad_names = names, bad_sequences = sequences,
    comments = comments
  )

  expect_error(dataset <- sequence_data_vector$new(sequences = bad_data))

  dataset <- sequence_data_vector$new()
  dataset$add_seqs(names, sequences, comments)

  df <- dataset$export()

  expect_equal(dataset$get_num_seqs(), 3)
  expect_equal(df$names, names)
  expect_equal(df$sequences, sequences)
  expect_equal(df$comments, comments)
  expect_equal(df$ends, c(49, 49, 36))
  expect_equal(df$starts, c(1, 1, 8))
  expect_equal(df$total_abundance, c(1, 1, 1))
  expect_equal(dataset$get_names()[[2]], "seq2")


  dataset$write_fasta_file(df$names, df$sequences,
    filename =
      "test_dataset_sequences.fasta"
  )

  # test initialization with file
  dataset <- sequence_data_vector$new(
    filename = "test_dataset_sequences.fasta",
    path = "."
  )
  expect_equal(dataset$get_num_seqs(), 3)

  # merge seqs without group data
  seqs_to_merge <- list(c("mySeq", "seq2", "seq3"))
  dataset$merge_seqs(seqs_to_merge)
  expect_equal(dataset$get_num_seqs(), 3)
  expect_equal(dataset$get_num_unique_seqs(), 1)

  dataset$clear()
  expect_equal(dataset$get_num_seqs(), 0)

  dataset$add_seqs(names, sequences, comments)
  expect_equal(dataset$get_num_seqs(), 3)

  remove_file("test_dataset_sequences.fasta")

  expect_equal(dataset$get_align_report(), data.table())
  expect_equal(dataset$get_contigs_report(), data.table())
})

test_that("test sequence_data_vector count info ", {
  sequences <- c(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    ".......GAGC-TTATCCGGATT-ATTG-GTT-AAA----",
    "...NNNGAGC-TTATCCGGATT-ATTG-GTT-AAA-GU-"
  )
  names <- c("seq1", "seq2", "seq3", "seq4")
  groups <- c("sample1", "sample2", "sample1", "sample2")

  dataset <- sequence_data_vector$new()
  dataset$add_seqs(names, sequences)
  expect_equal(dataset$get_num_seqs(), 4)
  expect_equal(dataset$get_group_totals(), 0)

  dataset$set_group_assignments(names, groups)

  expect_equal(dataset$get_groups(), c("sample1", "sample2"))
  expect_equal(dataset$get_names(group = "sample1"), c("seq1", "seq3"))
  expect_equal(dataset$get_names(group = "sample2"), c("seq2", "seq4"))
  expect_equal(dataset$get_names(group = "invalid_sample"), c(
    "seq1", "seq2",
    "seq3", "seq4"
  ))
  sample1_seqs <- c(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    ".......GAGC-TTATCCGGATT-ATTG-GTT-AAA----"
  )
  sample2_seqs <- c(
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    "...NNNGAGC-TTATCCGGATT-ATTG-GTT-AAA-GU-"
  )

  expect_equal(dataset$get_seqs(group = "sample1"), sample1_seqs)
  expect_equal(dataset$get_seqs(group = "sample2"), sample2_seqs)
  expect_equal(dataset$get_seqs(group = "invalid_sample"), sequences)
  expect_equal(dataset$get_seqs_abunds(bysample = TRUE)[[3]], c(1, 0))
  expect_equal(dataset$get_seqs_abunds(), c(1, 1, 1, 1))
  expect_equal(dataset$get_groups("seq1"), c("sample1"))
  expect_equal(dataset$get_num_seqs("sample1"), 2)
  expect_equal(dataset$get_num_seqs("sampleNotInDataset"), 0)
  expect_equal(dataset$get_group_totals(), c(2, 2))
  expect_equal(dataset$get_num_seqs(), 4)
  expect_equal(dataset$get_num_groups(), 2)

  df <- dataset$export()

  expect_equal(df$names, names)
  expect_equal(df$sequences, sequences)
  expect_equal(df$comments, c("", "", "", ""))
  expect_equal(df$ends, c(49, 49, 36, 38))
  expect_equal(df$starts, c(1, 1, 8, 4))
  expect_equal(df$sample1, c(1, 0, 1, 0))
  expect_equal(df$sample2, c(0, 1, 0, 1))

  dataset$clear()
  dataset$add_seqs(names, sequences)
  bad_groups <- c("sample1", "sample1", "sample2")
  bad_names <- c("nameNotPreset", "seq2", "seq3", "seq4")

  expect_error(dataset$set_group_assignments(names, bad_groups))
  expect_error(dataset$set_group_assignments(bad_names, groups))
  expect_error(dataset$set_abundance_assignments(bad_names, groups))
  expect_error(dataset$set_abundance_assignments(names, bad_groups))

  # adding seqs with abundance of 1
  dataset$clear()
  dataset$add_seqs(names, sequences)
  expect_equal(dataset$get_abund("seq2"), 1)
  expect_equal(dataset$get_abund("nameNotPreset"), 0)
  expect_equal(dataset$get_abunds("nameNotPreset"), c())
  expect_equal(dataset$get_groups("nameNotPreset"), c())


  abunds <- list(c(1150), c(115), c(50), c(4))
  dataset$clear()
  dataset$add_seqs(names, sequences)
  expect_equal(dataset$get_num_seqs(), 4)
  dataset$set_abundances(names, abunds)
  expect_equal(dataset$get_abund("seq2"), 115)
  expect_equal(dataset$get_num_seqs(), 1319)

  # dataset with no groups
  input <- c(
    "Representative_Sequence     total",
    "seq1	1150", "seq2	115", "seq3	50", "seq4	4"
  )

  create_dummy_file("test_no_groups.count_table", input)

  dataset <- sequence_data_vector$new()
  dataset$add_seqs(names, sequences)

  dataset$set_group_assignments(filename = "test_no_groups.count_table")

  remove_file("test_no_groups.count_table")

  expect_equal(dataset$get_num_seqs(), 1319)
  expect_equal(dataset$get_num_groups(), 0)
  dataset$write_count_file(
    "test_write_no_group.count_table"
  )

  dataset$clear()
  dataset$add_seqs(names, sequences)
  dataset$set_group_assignments(
    filename =
      "test_write_no_group.count_table"
  )
  expect_equal(dataset$get_num_seqs(), 1319)
  expect_equal(dataset$get_num_groups(), 0)

  remove_file("test_write_no_group.count_table")

  input <- c(
    "Representative_Sequence     total   sample2	sample3	sample4",
    "seq1	1150	250	400	500", "seq2	115	25	40	50",
    "seq3	50	25	25	0", "seq4	4	0	0	4"
  )

  create_dummy_file("test_small_abund.count_table", input)

  dataset$clear()
  dataset$add_seqs(names, sequences)
  dataset$set_group_assignments(
    filename =
      "test_small_abund.count_table"
  )
  expect_equal(dataset$get_num_seqs(), 1319)
  expect_equal(dataset$get_num_groups(), 3)

  remove_file("test_small_abund.count_table")

  dataset$write_count_file(
    "test_write.count_table"
  )

  dataset$clear()
  dataset$add_seqs(names, sequences)
  dataset$set_group_assignments(
    filename =
      "test_write.count_table"
  )
  expect_equal(dataset$get_num_seqs(), 1319)
  expect_equal(dataset$get_num_groups(), 3)

  remove_file("test_write.count_table")

  input <- c(
    "#Compressed Format: ...",
    "#1,F3D0	2,F3D1	3,F3D141	4,F3D142	5,F3D143",
    "Representative_Sequence	total	F3D0	F3D1	F3D141	F3D142	F3D143",
    "seq1	211	1,1 2,10 5,200",
    "seq2	100	4,100",
    "seq3	50	4,25 5,25",
    "seq4	115	1,100 3,10 5,5",
    "seq5	25	2,25",
    "seq6	1	3,1"
  )

  create_dummy_file("test_compressed.count_table", input)

  sequences <- c(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    ".......GAGC-TTATCCGGATT-ATTG-GTT-AAA----",
    "...NNNGAGC-TTATCCGGATT-ATTG-GTT-AAA-GU-",
    "ATGCCGT--CGAAT--CCGGTTAAA---TGT-GA...",
    "...NNTGCAA-ATGCCGT--CGAAT--CCGGTTAAA---TGT-GA..."
  )
  names <- c("seq1", "seq2", "seq3", "seq4", "seq5", "seq6")

  dataset$clear()
  dataset$add_seqs(names, sequences)
  dataset$set_group_assignments(
    filename =
      "test_compressed.count_table"
  )
  remove_file("test_compressed.count_table")

  expect_equal(dataset$get_num_seqs(), 502)
  expect_equal(dataset$get_num_groups(), 5)

  groups_totals <- dataset$get_group_totals()

  # 101  35  11 125 230
  expect_equal(groups_totals[1], 101)
  expect_equal(groups_totals[2], 35)
  expect_equal(groups_totals[3], 11)
  expect_equal(groups_totals[4], 125)
  expect_equal(groups_totals[5], 230)
  expect_equal(dataset$get_groups()[3], "F3D141")
  expect_equal(dataset$get_abund("seq4", "F3D141"), 10)
  expect_equal(dataset$get_abund("seq4", "F3D1"), 0)

  # remove all seqs from F3D141 to make sure group is removed as well
  names_to_remove <- c("seq4", "seq6")
  trash_codes <- c("testRemoveGroup_F3D141,", "testRemoveGroup_F3D141,")

  dataset$remove_seqs(names_to_remove, trash_codes)
  expect_equal(dataset$get_num_seqs(), 386)
  # 4 groups instead of 5
  expect_equal(dataset$get_num_groups(), 4)

  accnos_summary <- dataset$get_accnos_report()
  expect_equal(accnos_summary$total_bad, 116)
  expect_equal(accnos_summary$unique_bad, 2)
  expect_equal(
    accnos_summary$summary$trash_codes[[1]],
    "testRemoveGroup_F3D141"
  )
  groups_totals <- dataset$get_group_totals()

  # 1  35   0 125 225
  expect_equal(groups_totals[1], 1)
  expect_equal(groups_totals[2], 35)
  expect_equal(groups_totals[3], 125)
  expect_equal(groups_totals[4], 225)
  expect_equal(dataset$get_groups()[3], "F3D142")

  names_to_remove <- c("seq1")
  trash_codes <- c("testRemoveSeq1,")
  # removing seq1 will also remove group F3D0
  dataset$remove_seqs(names_to_remove, trash_codes)
  expect_equal(dataset$get_num_seqs(), 175)
  expect_equal(dataset$get_num_groups(), 3)

  accnos_summary <- dataset$get_accnos_report()

  expect_equal(accnos_summary$total_bad, 327)
  expect_equal(accnos_summary$unique_bad, 3)
  expect_equal(
    accnos_summary$summary$trash_codes[[1]],
    "testRemoveSeq1"
  )

  groups_totals <- dataset$get_group_totals()

  # 25 125  25
  expect_equal(groups_totals[1], 25)
  expect_equal(groups_totals[2], 125)
  expect_equal(groups_totals[3], 25)
  expect_equal(dataset$get_groups()[3], "F3D143")

  dataset$reinstate_seqs(c("testRemoveSeq1"))

  expect_equal(dataset$get_num_seqs(), 386)
  expect_equal(dataset$get_num_groups(), 4)

  accnos_summary <- dataset$get_accnos_report()
  expect_equal(accnos_summary$total_bad, 116)
  expect_equal(accnos_summary$unique_bad, 2)
  expect_equal(
    accnos_summary$summary$trash_codes[[1]],
    "testRemoveGroup_F3D141"
  )
  groups_totals <- dataset$get_group_totals()

  # 1  35   0 125 225
  expect_equal(groups_totals[1], 1)
  expect_equal(groups_totals[2], 35)
  expect_equal(groups_totals[3], 125)
  expect_equal(groups_totals[4], 225)
  expect_equal(dataset$get_num_seqs(), 386)
  expect_equal(dataset$get_groups()[3], "F3D142")

  # set seq1 and seq2 abundances
  expect_equal(dataset$get_abunds("seq1"), c(1L, 10L, 0L, 200L))
  expect_equal(dataset$get_abunds("seq2"), c(0L, 0L, 100L, 0L))

  # before setting abundance seq1 -> 1,10,0,200. seq2 -> 0,0,100,0

  dataset$set_abundances(c("seq1", "seq2"), list(c(0, 0, 30, 0), c(3, 4, 5, 6)))
  expect_equal(dataset$get_abunds("seq1"), c(0L, 0L, 30L, 0L))
  expect_equal(dataset$get_abunds("seq2"), c(3L, 4L, 5L, 6L))
  expect_equal(dataset$get_num_seqs(), 123)
  expect_equal(dataset$get_num_groups(), 4)

  groups_totals <- dataset$get_group_totals()

  # groups before set_abundances: 1  35   0 125 225
  # groups after set_abundances:  3  29.  0 65 31
  expect_equal(groups_totals[1], 3)
  expect_equal(groups_totals[2], 29)
  expect_equal(groups_totals[3], 60)
  expect_equal(groups_totals[4], 31)

  count_table <- dataset$get_count_table()
  expect_equal(count_table$F3D142[1], 30)
  expect_equal(count_table$F3D142[3], 25)

  seqs <- c("ATGCTGGG", "ATGGGGGG", "ATGATTTAG", "ATGGCCTCTA")
  dataset$set_seqs(dataset$get_names(), seqs)
  expect_equal(dataset$get_seqs(), seqs)

  expect_equal(dataset$get_num_unique_seqs(), 4)
  expect_equal(dataset$get_num_seqs(), 123)

  n <- c("seq1", "seq2", "seq3", "seq5")
  expect_equal(dataset$get_names(), n)

  seqs_to_merge <- list(c("seq1", "seq2", "seq3"), c("seq5"))
  dataset$merge_seqs(seqs_to_merge)

  expect_equal(dataset$get_num_unique_seqs(), 2)
  expect_equal(dataset$get_num_seqs(), 123)
  expect_equal(groups_totals[1], 3)
  expect_equal(groups_totals[2], 29)
  expect_equal(groups_totals[3], 60)
  expect_equal(groups_totals[4], 31)

  dataset$set_abundances(c("seq1", "seq5"), list(c(0, 0, 30, 0), c(0, 0, 0, 0)))
  expect_equal(dataset$get_num_unique_seqs(), 1)
  expect_equal(dataset$get_num_seqs(), 30)
  expect_equal(dataset$get_num_groups(), 1)
  groups_totals <- dataset$get_group_totals()
  expect_equal(groups_totals[1], 30)
})
