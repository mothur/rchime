# test sequence_dataset R6 Class

test_that("test sequence_dataset R6 class", {
  expect_error(dataset <- sequence_data_vector$new(sequences = "seq1"))

  dataset <- sequence_data_vector$new()
  df <- dataset$read_fasta_file(filename = rchime_example("test.fasta"))

  expect_equal(nrow(df), 323)
  expect_equal(df$names[[5]], "M00967_43_000000000-A3JHG_1_1111_23613_8947")
  expect_equal(df$starts[[5]], 1)
  expect_equal(df$ends[[5]], 278)
  expect_equal(df$polymers[[5]], 4)

  dataset$write_fasta_file(df$names, df$sequences,
    filename =
      "test_dataset_sequences.fasta"
  )

  # test initialization with file
  dataset <- sequence_data_vector$new(
    filename = "test_dataset_sequences.fasta",
    path = "."
  )
  remove_file("test_dataset_sequences.fasta")

  expect_equal(dataset$get_num_seqs(), 323)
  expect_equal(
    dataset$get_names()[5],
    "M00967_43_000000000-A3JHG_1_1111_23613_8947"
  )

  dataset$set_group_assignments(filename = rchime_example("test.count_table"))

  seqs <- dataset$get_seqs_table()
  expect_equal(seqs$names[[5]], "M00967_43_000000000-A3JHG_1_1111_23613_8947")
  expect_equal(substr(seqs$sequences[[5]], 1, 20), "TACGTAGGGG-GCAAGC-GT")
  expect_equal(dataset$get_group_totals()[1], 102970)

  dataset$clear()
  expect_equal(dataset$get_num_seqs(), 0)

  sequences <- c(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGT",
    ".......GAGC-TTATCCGGATT-ATTG-GTT-AAA----"
  )
  names <- c("mySeq", "seq2", "seq3")

  dataset$add_seqs(names, sequences)
  expect_equal(dataset$get_num_seqs(), 3)
  expect_equal(dataset$get_group_totals(), 0)
})
