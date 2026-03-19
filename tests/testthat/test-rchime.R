# test rchime function

test_that("test rchime by reference ", {

    fasta_data <- strollur::read_fasta(rchime_example("miseq.ng.fasta"))
    ref_data <- strollur::read_fasta(rchime_example("silva.v4.ng.fasta.gz"))

    data <- strollur::new_dataset("rchime reference example")
    strollur::add(data, table = fasta_data, type = "sequences")

    reference <- new_dataset("Silva V4 Region")
    strollur::add(reference, table = ref_data, type = "sequences")

    chimera_report <- rchime(data, reference = reference)

    # checks to make sure the correct things are created
    expect_equal(length(chimera_report), 2)
    expect_equal(length(chimera_report$accnos), 1676)
    expect_equal(nrow(chimera_report$rchime_report), 6084)
})

test_that("test rchime denovo - single sample ", {
    # no need to check dereplicate since one sample

    fasta_data <- strollur::read_fasta(rchime_example("miseq.ng.fasta"))
    count_file <- rchime_example("single_sample.count_table")
    abundance_data <- strollur::read_mothur_count(count_file)

    data <- strollur::new_dataset("rchime denovo example")
    strollur::add(data, table = fasta_data, type = "sequences")
    strollur::assign(data, table = abundance_data, type = "sequence_abundance")

    chimera_report <- rchime(data)

    # checks to make sure the correct things are created
    expect_equal(length(chimera_report), 2)
    expect_equal(length(chimera_report$accnos), 3719)
    expect_equal(nrow(chimera_report$rchime_report), 6084)
})

test_that("test rchime denovo - mulitple samples ", {

    fasta_data <- strollur::read_fasta(rchime_example("miseq.ng.fasta"))
    count_file <- rchime_example("miseq.count_table")
    abundance_data <- strollur::read_mothur_count(count_file)

    data <- strollur::new_dataset("rchime denovo example")
    strollur::add(data, table = fasta_data, type = "sequences")
    strollur::assign(data, table = abundance_data, type = "sequence_abundance")

    chimera_report <- rchime(data, dereplicate = TRUE, silent = TRUE)

    # checks to make sure the correct things are created
    expect_equal(length(chimera_report), 3)
    expect_equal(length(chimera_report$accnos), 3588)
    expect_equal(nrow(chimera_report$rchime_report), 6084)
})
