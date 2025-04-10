# test chimera_uchime function

test_that("test chimera_uchime denovo ", {
  names <- c(
    "chimera1", "non_chimera2ParentB", "chimera1ParentB",
    "non_chimera2", "chimera1ParentA", "non_chimera2ParentA"
  )

  chimera1_seq <- paste0("TACGG-AGGAT-GCGA-G-C-G-TTATCCGG-ATTTATTGGGTTT-AAA-GG",
    "GT-GCG-CAGGC-GG-ACTCTC-AA-GTCA-GCGGTCAAA-TCGC-GG-GG-CT",
    "C-AA-CC-CC-G-TT-CCG-CCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAA",
    "G-TAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAG-AT-ATCACGCAG-A",
    "ACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGA-GGCACGAA",
    "AGCGTGGGTATCGAACAGG",
    sep = ""
  )

  non_chimera2_pb <- paste0("TACGT-AGGGG-GCAA-G-C-G-TTATCCGG-ATTTACTGGGTGT-AAA",
    "-GGGA-GCG-TAGAC-GG-CAGCGC-AA-GTCT-GGAGTGAAA-TGCC-GG",
    "-GG-CCC-AA-CC-CC-G-GAACTG-CTTTGGAAACTGTGCAGCTCGAGT",
    "GCAGGAGAGG-TAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAG-AT",
    "-ATTAGGAGG-AACACCAGTGGCGAAGGCGGCTTACTGGACTGTAACTGA",
    "CGTTGA-GGCTCGAAAGCGTGGGGAGCAAACAGG",
    sep = ""
  )

  chimera1_pb_seq <- paste0("TACGG-AGGAT-GCGA-G-C-G-TTATCCGG-ATTTATTGGGTTT-AAA",
    "-GGGT-GCG-TAGGC-GG-GCTGTT-AA-GTCA-GCGGTCAAA-TGTC-GG-G",
    "G-CTC-AA-CC-CC-G-GC-CTG-CCGTTGAAACTGGCGGCCTCGAGTGGG",
    "CGAGAAG-TATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAG-AT-ATC",
    "ACGCAG-AACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTG",
    "A-GGCACGAAAGCGTGGGTATCGAACAGG",
    sep = ""
  )

  non_chimera2_seq <- paste0("TACGT-AGGGG-GCAA-G-C-G-TTATCCGG-ATTTACTGGGTGT-AA",
    "A-GGGA-GCG-TAGGC-GG-CCATGC-AA-GTCA-GAAGTGAAA-ACCC-",
    "GG-GG-CTC-AA-C--CC-TGGGAGTG-CTTTTGAAACTGTGCGGCTAG",
    "AGTGTCGGAGGGG-TAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTA",
    "G-AT-ATTAGGAGG-AACACCAGTGGCGAAGGCGGCTTACTGGACGATG",
    "ACTGACGCTGA-GGCTCGAAAGCGTGGGGAGCAAACAGG",
    sep = ""
  )

  chimera1_pa_seq <- paste0("TACGG-AGGAT-GCGA-G-C-G-TTATCCGG-ATTTATTGGGTTT-AAA",
    "-GGGT-GCG-CAGGC-GG-ACTCTC-AA-GTCA-GCGGTCAAA-TCGC-GG-G",
    "G-CTC-AA-CC-CC-G-TT-CCG-CCGTTGAAACTGGGAGCCTTGAGTGCG",
    "CGAGAAG-TAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAG-AT-ATC",
    "ACGCAG-AACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTC",
    "A-TGCACGAAAGCGTGGGTATCGAACAGG",
    sep = ""
  )

  non_chimera2_pa_seq <- paste0("TACGT-AGGGG-GCAA-G-C-G-TTATCCGG-ATTTACTGGGTGT",
    "-AAA-GGGA-GCG-TAGAC-GG-TTTTGC-AA-GTCT-GATGTGAAA-GCCC-GG",
    "-GG-CTC-AA-CC-CC-G-GGACTG-CATTGGAAACTGTATGACTGGAGT",
    "GCAGGAGAGG-TAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAG-AT",
    "-ATTAGGAGG-AACACCAGTGGCGAAGGCGGCTTACTGGACTGTAACTGA",
    "CGTTGA-GGCTCGAAAGCGTGGGGAGCAAACAGG",
    sep = ""
  )


  seqs <- c(
    chimera1_seq, non_chimera2_pb, chimera1_pb_seq,
    non_chimera2_seq, chimera1_pa_seq, non_chimera2_pa_seq
  )

  dataset <- sequence_data_table$new()
  dataset$add_seqs(names, seqs)
  dataset$set_group_assignments(
    filename =
      rchime_example("test.chimera.count_table")
  )

  expect_equal(dataset$get_num_unique_seqs(), 6)
  expect_equal(dataset$get_num_seqs(), 295696)

  results <- chimera_uchime(dataset, dereplicate = TRUE, chimealns = TRUE)

  expect_equal(dataset$get_num_unique_seqs(), 5)
  expect_equal(dataset$get_num_seqs(), 290923)

  # results$accnos results
  expect_equal(results$accnos[[1]], "chimera1")

  # results$uchimeout - chimera1 results
  expect_equal(results$uchimeout$Query[5], "chimera1/ab=349/")
  expect_equal(results$uchimeout$A[5], "chimera1ParentA/ab=7912/")
  expect_equal(results$uchimeout$B[5], "chimera1ParentB/ab=5770/")
  expect_equal(results$uchimeout$IdQM[5], 100.0)
  expect_equal(round(results$uchimeout$IdQA[5], 2), 97.98)
  expect_equal(round(results$uchimeout$IdQB[5], 2), 94.35)
  expect_equal(round(results$uchimeout$IdAB[5], 2), 92.34)
  expect_equal(round(results$uchimeout$IdQT[5], 2), 98.02)
  expect_equal(round(results$uchimeout$Div[5], 2), 2.06)
  expect_equal(round(results$uchimeout$Score[5], 2), 0.56)
  expect_equal(results$uchimeout$LY[5], 14)
  expect_equal(results$uchimeout$LN[5], 0)
  expect_equal(results$uchimeout$LA[5], 0)
  expect_equal(results$uchimeout$RY[5], 5)
  expect_equal(results$uchimeout$RN[5], 0)
  expect_equal(results$uchimeout$RA[5], 0)
  expect_equal(results$uchimeout$Y[5], "Y")

  # results$uchimeout - non_chimera2 results
  expect_equal(results$uchimeout$Query[6], "non_chimera2/ab=328/")
  expect_equal(results$uchimeout$A[6], "non_chimera2ParentA/ab=1013/")
  expect_equal(results$uchimeout$B[6], "non_chimera2ParentB/ab=838/")
  expect_equal(round(results$uchimeout$IdQM[6], 2), 91.97)
  expect_equal(round(results$uchimeout$IdQA[6], 2), 91.16)
  expect_equal(round(results$uchimeout$IdQB[6], 2), 90.76)
  expect_equal(round(results$uchimeout$IdAB[6], 2), 93.17)
  expect_equal(round(results$uchimeout$IdQT[6], 2), 91.30)
  expect_equal(round(results$uchimeout$Div[6], 2), 0.88)
  expect_equal(round(results$uchimeout$Score[6], 2), 0.02)
  expect_equal(results$uchimeout$LY[6], 5)
  expect_equal(results$uchimeout$LN[6], 2)
  expect_equal(results$uchimeout$LA[6], 6)
  expect_equal(results$uchimeout$RY[6], 4)
  expect_equal(results$uchimeout$RN[6], 2)
  expect_equal(results$uchimeout$RA[6], 10)
  expect_equal(results$uchimeout$Y[6], "N")

  # results$uchimealns - chimera1 results
  expect_equal(results$uchimealns[[1]][1], "Query ( 252nt ) chimera1")
  expect_equal(results$uchimealns[[1]][2], "ParentA ( 252nt ) chimera1ParentA")
  expect_equal(results$uchimealns[[1]][3], "ParentB ( 252nt ) chimera1ParentB")

  model <- paste0("Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
    "BBBBBBBBBBBBBBBBBBBBBBBB",
    sep = ""
  )
})

test_that("test chimera_uchime reference ", {
  dataset <- sequence_data_table$new(
    filename = rchime_example("test.fasta")
  )

  reference <- sequence_data_table$new(
    filename = rchime_example("reference.fasta")
  )
})
