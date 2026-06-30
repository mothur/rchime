#' create internal strollur objects for tests and examples
#'
#' This function generates and saves the .rds strollur objects to inst/extdata.
#' @return Invisible NULL
generate_strollur_objects <- function() {

    output_dir <- "/Users/swestcot/Desktop/rchime/inst/extdata"

    # strollur_reference.rds
    filename <- file.path(output_dir, "strollur_reference.rds")
    data <- strollur::new_dataset("silva gold")
    strollur::add(data, silva_gold())
    strollur::save_dataset(data, filename)

    # strollur_multi_sample_small.rds - 500 sequences
    filename <- file.path(output_dir, "strollur_multi_sample_small.rds")
    table <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))
    data <- strollur::new_dataset("multi_sample - 500 sequences")
    strollur::add(data, table = unique(table[, c("sequence_name", "sequence")]),
                  type = "sequence")
    strollur::assign(data, table, type = "sequence_abundance")
    strollur::save_dataset(data, filename)

    # strollur_miseq_tiny.rds - 100 sequences - no abunds
    filename <- file.path(output_dir, "strollur_miseq_tiny.rds")
    table <- readRDS(rchime_example("miseq_fasta.rds"))
    data <- strollur::new_dataset("miseq - 100 sequences")
    strollur::add(data, table = table[1:100, ])
    strollur::save_dataset(data, filename)

    # strollur_single_sample_tiny.rds
    filename <- file.path(output_dir, "strollur_single_sample_tiny.rds")
    table <- readRDS(rchime_example("miseq_data_frame.rds"))
    data <- strollur::new_dataset("single_sample - 100 sequences")
    strollur::add(data, table = table[1:100, ])
    strollur::assign(data, table[1:100, ], type = "sequence_abundance")
    strollur::save_dataset(data, filename)

    invisible(NULL)
}
