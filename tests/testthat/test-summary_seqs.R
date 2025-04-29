# test summary_seqs function

test_that("test summary_seqs summarize fasta / count ", {
  # no valid inputs
  expect_error(summary_seqs())

  samples <- c(
    "sample1", "sample1", "sample1",
    "sample2", "sample3", "sample3"
  )
  sequences <- c(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    ".......GAGC-TTATCCGGATT-ATTG-GTT-AAA----",
    "NACGG-GGATGCG---TATCCGGATTTAT---TTTAAAGGGTGCGT",
    ".....GGATGCGAG-GTTATCC--ATTTATTGGGTTTA..........",
    ".TTTTTTTTTTTTTGAGC-TTATCCGGATT-ATTG-GTT-AAA----"
  )
  names <- c("seq1", "seq2", "seq3", "seq4", "seq5", "seq6")


  input <- c(
    "Representative_Sequence     total   sample2	sample3	sample4",
    "seq1	1150	250	400	500", "seq2	115	25	40	50",
    "seq3	50	25	25	0", "seq4	4	0	0	4",
    "seq5	1	1	0	0", "seq6	100	50	0	50"
  )

  create_dummy_file("test.count_table", input)

  dataset <- sequence_data_vector$new()

  # can't summarize an empty dataset
  expect_error(summary_seqs(dataset = dataset))

  dataset$add_seqs(names = names, sequences = sequences)

  # no count data
  results <- summary_seqs(dataset = dataset, silent = TRUE)

  results$fasta_summary <- setDF(results$fasta_summary)

  # 25% numseqs
  expect_equal(results$fasta_summary[3, 7], "2")
  # 75% numseqs
  expect_equal(results$fasta_summary[5, 7], "5")

  # add count data
  dataset2 <- sequence_data_vector$new()
  dataset2$add_seqs(names = names, sequences = sequences)
  dataset2$set_group_assignments(filename = "test.count_table")

  remove_file("test.count_table")

  all_results <- summary_seqs(dataset = dataset2)

  results <- setDF(all_results$fasta_summary)

  # minimum start
  expect_equal(results[1, 1], 1)
  # maximum start
  expect_equal(results[7, 1], 8)

  # minimum end
  expect_equal(results[1, 2], 36)
  # maximum end
  expect_equal(results[7, 2], 49)

  # minimum length
  expect_equal(results[1, 3], 25)
  # maximum length
  expect_equal(results[7, 3], 49)

  # 2.5% ambigs
  expect_equal(results[2, 4], 0)
  # 97.5% ambigs
  expect_equal(results[6, 4], 1)

  # 2.5% polymers
  expect_equal(results[2, 5], 3)
  # 97.5% polymers
  expect_equal(results[6, 5], 13)

  # 25% numseqs
  expect_equal(results[3, 7], "356")
  # 75% numseqs
  expect_equal(results[5, 7], "1066")

  # Means start, end, length, ambigs, polymer, num_ns
  expect_equal(results[8, 1], 1)
  expect_equal(results[8, 2], 48)
  expect_equal(results[8, 3], 47)
  expect_equal(results[8, 4], 0)
  expect_equal(results[8, 5], 3)
  expect_equal(results[8, 6], 0)
})
