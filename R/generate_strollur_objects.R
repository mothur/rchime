#' create internal strollur objects for tests and examples
#'
#' This function generates and saves the .rds strollur objects to inst/extdata.
#' @return Invisible NULL
#' @noRd
#' @keywords internal
generate_strollur_objects <- function() {
  output_dir <- "/Users/swestcot/Desktop/rchime/inst/extdata"

  # strollur_miseq_tiny.rds
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  filename <- file.path(output_dir, "strollur_miseq_tiny.rds")
  data <- strollur::new_dataset("100 sequences")
  strollur::xdev_add_sequences(data, table = fasta_data[1:100, ])
  strollur::save_dataset(data, filename)

  # strollur_reference.rds
  reference_data <- silva_gold()
  filename <- file.path(output_dir, "strollur_reference.rds")
  data <- strollur::new_dataset("reference sequences")
  strollur::xdev_add_sequences(data, table = reference_data)
  strollur::save_dataset(data, filename)

  # "strollur_multi_sample_small.rds"
  table <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))
  filename <- file.path(output_dir, "strollur_multi_sample_small.rds")
  data <- strollur::new_dataset("miseq - 500 sequences")
  strollur::xdev_add_sequences(data, table = unique(table[, c(
    "sequence_name",
    "sequence"
  )]))
  strollur::assign(data, table = table, type = "sequence_abundance")
  strollur::save_dataset(data, filename)

  # strollur_single_sample.rds
  table <- readRDS(rchime_example("miseq_data_frame.rds"))
  filename <- file.path(output_dir, "strollur_single_sample.rds")
  data <- strollur::new_dataset("single sample")
  strollur::xdev_add_sequences(data, table = table)
  strollur::assign(data, table = table, type = "sequence_abundance")
  strollur::save_dataset(data, filename)

  # strollur_multi_sample_tiny.rds
  table <- readRDS(rchime_example("miseq_data_frame.rds"))
  filename <- file.path(output_dir, "strollur_multi_sample_tiny.rds")
  data <- strollur::new_dataset("miseq - 500 sequences")
  strollur::xdev_add_sequences(data, table = table)
  strollur::assign(data, table = table, type = "sequence_abundance")
  strollur::save_dataset(data, filename)

  invisible(NULL)
}
