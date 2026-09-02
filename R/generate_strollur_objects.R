#' @return strollur::strollur object with 100 sequences
#' @noRd
#' @keywords internal
strollur_miseq_tiny <- function() {
  fasta_data <- readRDS(rchime_example("miseq_fasta.rds"))
  data <- strollur::new_dataset("100 sequences")
  strollur::xdev_add_sequences(data, table = fasta_data[1:100, ])
}

#' @return strollur::strollur object with silva gold sequences
#' @noRd
#' @keywords internal
strollur_reference <- function() {
  data <- strollur::new_dataset("silva gold reference sequences")
  strollur::xdev_add_sequences(data, table = silva_gold())
}

#' @title Create a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur object}
#'   object with 500 sequences
#' @description
#' Create a \href{https://mothur.org/strollur/reference/strollur.html}{strollur
#'   object} object with 500 sequences
#' @examples
#' rchime::strollur_multi_sample_small()
#' @return \href{https://mothur.org/strollur/reference/strollur.html}{strollur
#'   object} with 500 sequences
#' @export
strollur_multi_sample_small <- function() {
  table <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))
  data <- strollur::new_dataset("miseq - 500 sequences")
  strollur::xdev_add_sequences(data, table = unique(table[, c(
    "sequence_name",
    "sequence"
  )]))
  strollur::assign(data, table = table, type = "sequence_abundance")
}

#' @return strollur::strollur object with single sample
#' @noRd
#' @keywords internal
strollur_single_sample <- function() {
  table <- readRDS(rchime_example("miseq_data_frame.rds"))
  data <- strollur::new_dataset("single sample")
  strollur::xdev_add_sequences(data, table = table)
  strollur::assign(data, table = table, type = "sequence_abundance")
}

#' @return strollur::strollur object with 500 sequences and multiple samples
#' @noRd
#' @keywords internal
strollur_multi_sample_tiny <- function() {
  table <- readRDS(rchime_example("miseq_data_frame.rds"))
  data <- strollur::new_dataset("miseq - 500 sequences")
  strollur::xdev_add_sequences(data, table = table)
  strollur::assign(data, table = table, type = "sequence_abundance")
}
