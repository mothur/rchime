# test rchime function

test_that("test uchime parameters ", {
  # minh < score == chimeric,
  # setting minh to 0.6 should cause chimera1 to be non_chimeric
  options <- rchime_options(
    minh = 0.6, abskew = 1.9, mindiv = 0.505,
    xn = 8.01, dn = 1.401, xa = 1.001,
    minchunk = 63, chunks = 5, idsmoothwindow = 33,
    maxp = 3,
    skipgaps = FALSE, skipgaps2 = FALSE
  )

  results <- rchime(
    fasta = rchime_example("test.chimera.fasta"),
    count = rchime_example("test_nogroups.count_table"),
    dereplicate = TRUE, silent = TRUE,
    rchime_options = options
  )

  # results$accnos results
  expect_equal(length(results$accnos), 0)

  # results$uchimeout - chimera1 results
  expect_equal(results$uchimeout$Query[5], "chimera1/ab=4773/")
  expect_equal(results$uchimeout$A[5], "chimera1ParentA/ab=150122/")
  expect_equal(results$uchimeout$B[5], "chimera1ParentB/ab=108939/")
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
  expect_equal(results$uchimeout$Y[5], "N")

  # results$uchimeout - non_chimera2 results
  expect_equal(results$uchimeout$Query[6], "non_chimera2/ab=1734/")
  expect_equal(results$uchimeout$A[6], "non_chimera2ParentA/ab=15750/")
  expect_equal(results$uchimeout$B[6], "non_chimera2ParentB/ab=14378/")
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
})

test_that("test rchime errors ", {
  # no inputs
  expect_error(rchime())

  dataset <- sequence_data_vector$new(rchime_example("test.chimera.fasta"))

  # too many inputs
  expect_error(rchime(
    dataset = dataset,
    fasta = rchime_example("test.chimera.fasta"),
    count = rchime_example("test_nogroups.count_table")
  ))

  # bad fasta file
  expect_error(rchime(fasta = "invalid_file.txt"))

  # unaligned dataset
  dataset$add_seqs(c("bad_seq"), c("ATGC-TTG"))
  expect_error(rchime(dataset))

  # empty dataset
  dataset$clear()
  expect_error(rchime(dataset))
})

test_that("test rchime denovo - no groups", {
  results <- rchime(
    fasta = rchime_example("test.chimera.fasta"),
    count = rchime_example("test_nogroups.count_table"),
    dereplicate = TRUE, chimealns = TRUE, silent = TRUE
  )

  expect_equal(results$dataset$get_num_unique_seqs(), 5)
  expect_equal(results$dataset$get_num_seqs(), 290923)

  # results$accnos results
  expect_equal(results$accnos[[1]], "chimera1")

  # results$uchimeout - chimera1 results
  expect_equal(results$uchimeout$Query[5], "chimera1/ab=4773/")
  expect_equal(results$uchimeout$A[5], "chimera1ParentA/ab=150122/")
  expect_equal(results$uchimeout$B[5], "chimera1ParentB/ab=108939/")
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
  expect_equal(results$uchimeout$Query[6], "non_chimera2/ab=1734/")
  expect_equal(results$uchimeout$A[6], "non_chimera2ParentA/ab=15750/")
  expect_equal(results$uchimeout$B[6], "non_chimera2ParentB/ab=14378/")
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

test_that("test rchime denovo - groups, dereplicate = TRUE / FALSE", {
  dataset <- sequence_data_vector$new(rchime_example("test.chimera.fasta"))
  dataset$set_group_assignments(
    filename =
      rchime_example("test.chimera.count_table")
  )

  expect_equal(dataset$get_num_unique_seqs(), 6)
  expect_equal(dataset$get_num_seqs(), 295696)

  results <- rchime(dataset, dereplicate = TRUE, chimealns = TRUE)

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

  dataset <- sequence_data_vector$new(rchime_example("test.chimera.fasta"))
  dataset$set_group_assignments(
    filename =
      rchime_example("test.chimera.count_table")
  )

  expect_equal(dataset$get_num_unique_seqs(), 6)
  expect_equal(dataset$get_num_seqs(), 295696)

  results <- rchime(dataset, dereplicate = FALSE, silent = TRUE)

  expect_equal(dataset$get_num_unique_seqs(), 5)
  expect_equal(dataset$get_num_seqs(), 290923)

  # results$accnos results
  expect_equal(results$accnos[[1]], "chimera1")
})


test_that("test rchime reference ", {
  names <- c("chimera1", "non_chimera2")

  chimera1_seq <- paste0("TACGG-AGGAT-GCGA-G-C-G-TTATCCGG-ATTTATTGGGTTT-AAA-GG",
    "GT-GCG-CAGGC-GG-ACTCTC-AA-GTCA-GCGGTCAAA-TCGC-GG-GG-CT",
    "C-AA-CC-CC-G-TT-CCG-CCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAA",
    "G-TAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAG-AT-ATCACGCAG-A",
    "ACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGA-GGCACGAA",
    "AGCGTGGGTATCGAACAGG",
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

  seqs <- c(chimera1_seq, non_chimera2_seq)

  dataset <- sequence_data_vector$new()
  dataset$add_seqs(names, seqs)

  expect_equal(dataset$get_num_unique_seqs(), 2)

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

  ref_seqs <- c(
    ">non_chimera2_pb", non_chimera2_pb,
    ">chimera1_pb_seq", chimera1_pb_seq,
    ">chimera1_pa_seq", chimera1_pa_seq,
    ">non_chimera2_pa_seq", non_chimera2_pa_seq
  )

  create_dummy_file("reference.fasta", ref_seqs)

  results <- rchime(dataset,
    reference = "reference.fasta",
    dereplicate = TRUE, chimealns = TRUE
  )

  remove_file("reference.fasta")

  expect_equal(dataset$get_num_unique_seqs(), 1)

  # results$accnos results
  expect_equal(results$accnos[[1]], "chimera1")

  # results$uchimeout - chimera1 results
  expect_equal(results$uchimeout$Query[1], "chimera1/ab=1/")
  expect_equal(results$uchimeout$A[1], "chimera1_pa_seq/ab=1/")
  expect_equal(results$uchimeout$B[1], "chimera1_pb_seq/ab=1/")
  expect_equal(results$uchimeout$IdQM[1], 100.0)
  expect_equal(round(results$uchimeout$IdQA[1], 2), 97.98)
  expect_equal(round(results$uchimeout$IdQB[1], 2), 94.35)
  expect_equal(round(results$uchimeout$IdAB[1], 2), 92.34)
  expect_equal(round(results$uchimeout$IdQT[1], 2), 98.02)
  expect_equal(round(results$uchimeout$Div[1], 2), 2.06)
  expect_equal(round(results$uchimeout$Score[1], 2), 0.56)
  expect_equal(results$uchimeout$LY[1], 14)
  expect_equal(results$uchimeout$LN[1], 0)
  expect_equal(results$uchimeout$LA[1], 0)
  expect_equal(results$uchimeout$RY[1], 5)
  expect_equal(results$uchimeout$RN[1], 0)
  expect_equal(results$uchimeout$RA[1], 0)
  expect_equal(results$uchimeout$Y[1], "Y")

  # results$uchimeout - non_chimera2 results
  expect_equal(results$uchimeout$Query[2], "non_chimera2/ab=1/")
  expect_equal(results$uchimeout$A[2], "non_chimera2_pa_seq/ab=1/")
  expect_equal(results$uchimeout$B[2], "non_chimera2_pb/ab=1/")
  expect_equal(round(results$uchimeout$IdQM[2], 2), 91.97)
  expect_equal(round(results$uchimeout$IdQA[2], 2), 91.16)
  expect_equal(round(results$uchimeout$IdQB[2], 2), 90.76)
  expect_equal(round(results$uchimeout$IdAB[2], 2), 93.17)
  expect_equal(round(results$uchimeout$IdQT[2], 2), 91.30)
  expect_equal(round(results$uchimeout$Div[2], 2), 0.88)
  expect_equal(round(results$uchimeout$Score[2], 2), 0.02)
  expect_equal(results$uchimeout$LY[2], 5)
  expect_equal(results$uchimeout$LN[2], 2)
  expect_equal(results$uchimeout$LA[2], 6)
  expect_equal(results$uchimeout$RY[2], 4)
  expect_equal(results$uchimeout$RN[2], 2)
  expect_equal(results$uchimeout$RA[2], 10)
  expect_equal(results$uchimeout$Y[2], "N")

  # results$uchimealns - chimera1 results
  expect_equal(results$uchimealns[[1]][1], "Query ( 252nt ) chimera1")
  expect_equal(results$uchimealns[[1]][2], "ParentA ( 252nt ) chimera1_pa_seq")
  expect_equal(results$uchimealns[[1]][3], "ParentB ( 252nt ) chimera1_pb_seq")

  model <- paste0("Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
    "BBBBBBBBBBBBBBBBBBBBBBBB",
    sep = ""
  )
})

test_that("test rchime compare mothur ", {

    mothur_accnos_deplicate_true <- c(
        "M00967_43_000000000-A3JHG_1_2106_11200_6323",
        "M00967_43_000000000-A3JHG_1_2109_12933_23855",
        "M00967_43_000000000-A3JHG_1_1113_13160_8353",
        "M00967_43_000000000-A3JHG_1_2103_18464_5156",
        "M00967_43_000000000-A3JHG_1_1104_19228_20175",
        "M00967_43_000000000-A3JHG_1_1111_9686_17749",
        "M00967_43_000000000-A3JHG_1_2102_16849_3037",
        "M00967_43_000000000-A3JHG_1_1113_6703_7749",
        "M00967_43_000000000-A3JHG_1_2101_9573_6498",
        "M00967_43_000000000-A3JHG_1_2104_21607_19035",
        "M00967_43_000000000-A3JHG_1_1111_11333_11888",
        "M00967_43_000000000-A3JHG_1_1108_16580_11543",
        "M00967_43_000000000-A3JHG_1_1109_20720_4988",
        "M00967_43_000000000-A3JHG_1_1111_25827_17738",
        "M00967_43_000000000-A3JHG_1_1106_27493_9675",
        "M00967_43_000000000-A3JHG_1_2113_18982_5442",
        "M00967_43_000000000-A3JHG_1_2113_10256_22364",
        "M00967_43_000000000-A3JHG_1_2106_26339_22698",
        "M00967_43_000000000-A3JHG_1_1112_21347_13848",
        "M00967_43_000000000-A3JHG_1_2104_20483_7812",
        "M00967_43_000000000-A3JHG_1_1111_14999_28567",
        "M00967_43_000000000-A3JHG_1_1114_26040_9558",
        "M00967_43_000000000-A3JHG_1_2112_5637_20821",
        "M00967_43_000000000-A3JHG_1_2107_15431_18507",
        "M00967_43_000000000-A3JHG_1_1106_11047_22258",
        "M00967_43_000000000-A3JHG_1_2113_22396_11112",
        "M00967_43_000000000-A3JHG_1_1113_12090_13717",
        "M00967_43_000000000-A3JHG_1_2113_28267_11188",
        "M00967_43_000000000-A3JHG_1_1105_10343_20641",
        "M00967_43_000000000-A3JHG_1_1101_20829_15688",
        "M00967_43_000000000-A3JHG_1_1105_11964_9315",
        "M00967_43_000000000-A3JHG_1_1113_24481_7796",
        "M00967_43_000000000-A3JHG_1_2104_15203_22079",
        "M00967_43_000000000-A3JHG_1_2104_10360_6292",
        "M00967_43_000000000-A3JHG_1_1114_12039_13282",
        "M00967_43_000000000-A3JHG_1_1114_17672_13068",
        "M00967_43_000000000-A3JHG_1_2101_16885_12036")

    mothur_accnos_deplicate_true <- sort(mothur_accnos_deplicate_true)

    results <- rchime(
        fasta = rchime_example("test.fasta"),
        count = rchime_example("test.count_table"),
        dereplicate = TRUE, chimealns = TRUE, silent = TRUE
    )

    expect_equal(sort(results$accnos), mothur_accnos_deplicate_true)

    mothur_accnos_deplicate_false <- c(
      "M00967_43_000000000-A3JHG_1_2106_26339_22698",
      "M00967_43_000000000-A3JHG_1_2112_5637_20821",
      "M00967_43_000000000-A3JHG_1_1106_11047_22258",
      "M00967_43_000000000-A3JHG_1_1114_12039_13282",
      "M00967_43_000000000-A3JHG_1_1114_17672_13068",
      "M00967_43_000000000-A3JHG_1_1105_11964_9315",
      "M00967_43_000000000-A3JHG_1_1101_20829_15688",
      "M00967_43_000000000-A3JHG_1_2101_16885_12036",
      "M00967_43_000000000-A3JHG_1_1112_21347_13848",
      "M00967_43_000000000-A3JHG_1_1101_13440_22261",
      "M00967_43_000000000-A3JHG_1_1111_14999_28567",
      "M00967_43_000000000-A3JHG_1_2114_14891_13095",
      "M00967_43_000000000-A3JHG_1_2107_15431_18507",
      "M00967_43_000000000-A3JHG_1_2113_10256_22364",
      "M00967_43_000000000-A3JHG_1_1114_26040_9558",
      "M00967_43_000000000-A3JHG_1_2112_20262_16676",
      "M00967_43_000000000-A3JHG_1_1106_27493_9675",
      "M00967_43_000000000-A3JHG_1_1111_8959_4660",
      "M00967_43_000000000-A3JHG_1_1113_12090_13717",
      "M00967_43_000000000-A3JHG_1_1112_15471_15261",
      "M00967_43_000000000-A3JHG_1_2104_4953_11390",
      "M00967_43_000000000-A3JHG_1_2104_21607_19035",
      "M00967_43_000000000-A3JHG_1_2104_10360_6292",
      "M00967_43_000000000-A3JHG_1_2113_14237_26547",
      "M00967_43_000000000-A3JHG_1_2104_20483_7812",
      "M00967_43_000000000-A3JHG_1_2113_18982_5442",
      "M00967_43_000000000-A3JHG_1_2102_16849_3037",
      "M00967_43_000000000-A3JHG_1_1105_21012_11572",
      "M00967_43_000000000-A3JHG_1_1111_25827_17738",
      "M00967_43_000000000-A3JHG_1_2107_7301_13458",
      "M00967_43_000000000-A3JHG_1_1109_20720_4988",
      "M00967_43_000000000-A3JHG_1_1110_17575_13842",
      "M00967_43_000000000-A3JHG_1_1108_16580_11543",
      "M00967_43_000000000-A3JHG_1_2113_28267_11188",
      "M00967_43_000000000-A3JHG_1_1111_11333_11888",
      "M00967_43_000000000-A3JHG_1_1105_10343_20641",
      "M00967_43_000000000-A3JHG_1_1112_16718_13888",
      "M00967_43_000000000-A3JHG_1_1113_6703_7749",
      "M00967_43_000000000-A3JHG_1_1101_19881_16112",
      "M00967_43_000000000-A3JHG_1_1111_26907_16090",
      "M00967_43_000000000-A3JHG_1_1114_10616_10969",
      "M00967_43_000000000-A3JHG_1_2114_22671_19859",
      "M00967_43_000000000-A3JHG_1_1105_24426_14439",
      "M00967_43_000000000-A3JHG_1_1111_9686_17749",
      "M00967_43_000000000-A3JHG_1_1113_24481_7796",
      "M00967_43_000000000-A3JHG_1_2113_22396_11112",
      "M00967_43_000000000-A3JHG_1_2109_12929_2238",
      "M00967_43_000000000-A3JHG_1_2104_15203_22079",
      "M00967_43_000000000-A3JHG_1_2110_19050_7165",
      "M00967_43_000000000-A3JHG_1_1113_3668_18339",
      "M00967_43_000000000-A3JHG_1_1114_16114_2126",
      "M00967_43_000000000-A3JHG_1_1104_19228_20175",
      "M00967_43_000000000-A3JHG_1_2103_18464_5156",
      "M00967_43_000000000-A3JHG_1_2101_9573_6498",
      "M00967_43_000000000-A3JHG_1_1113_13160_8353",
      "M00967_43_000000000-A3JHG_1_2111_17495_8499",
      "M00967_43_000000000-A3JHG_1_1111_12712_25853",
      "M00967_43_000000000-A3JHG_1_2109_12933_23855",
      "M00967_43_000000000-A3JHG_1_1114_24213_9688",
      "M00967_43_000000000-A3JHG_1_1108_23982_23453",
      "M00967_43_000000000-A3JHG_1_1104_21399_13763",
      "M00967_43_000000000-A3JHG_1_1102_24524_7196",
      "M00967_43_000000000-A3JHG_1_2109_16377_13199",
      "M00967_43_000000000-A3JHG_1_2111_10300_3779",
      "M00967_43_000000000-A3JHG_1_2114_8808_8472",
      "M00967_43_000000000-A3JHG_1_1107_10738_15216",
      "M00967_43_000000000-A3JHG_1_2106_11200_6323",
      "M00967_43_000000000-A3JHG_1_1104_22125_9596",
      "M00967_43_000000000-A3JHG_1_2114_21028_11825",
      "M00967_43_000000000-A3JHG_1_1110_18424_5907",
      "M00967_43_000000000-A3JHG_1_1104_14134_22528",
      "M00967_43_000000000-A3JHG_1_2112_28591_10042",
      "M00967_43_000000000-A3JHG_1_1104_4485_14998",
      "M00967_43_000000000-A3JHG_1_1105_19777_15422")

      mothur_accnos_deplicate_false <- sort(mothur_accnos_deplicate_false)

      results <- rchime(
          fasta = rchime_example("test.fasta"),
          count = rchime_example("test.count_table"),
          dereplicate = FALSE, silent = TRUE
      )

      expect_equal(sort(results$accnos), mothur_accnos_deplicate_false)

})
