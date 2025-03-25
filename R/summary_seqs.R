#' summary_seqs options
#'
#' @param dataset R6 'sequence_dataset' object containing sequence data
#' @param processors Integer, number of cores to use. Default = all available
#' @param silent Boolean, suppress progress updates and console output

#' @return list of data.tables containing the summary data
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @import cli
#' @importFrom parallelly, availableCores
#'
#' @examples
#'
#' # Let's create a dataset to summarize
#'
#' dataset <- sequence_data_table$new(
#'   filename = mothur2_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = mothur2_example("test.count_table")
#' )
#'
#' # To summarize the dataset, run the following:
#'
#' summary <- summary_seqs(dataset = dataset$dataset)
#'
#' @export
summary_seqs <- function(dataset = NULL,
                         processors = parallelly::availableCores(),
                         silent = FALSE) {
  if (is.null(dataset)) {
    message <- "[ERROR]: dataset is NULL. summary_seqs requires a
       'sequence_dataset' object."
    cli::cli_abort(message)
  }

  results_row_names <- c(
    "Minimum:", "2.5%-tile:", "25%-tile:",
    "Median:   ", "75%-tile:", "97.5%-tile:",
    "Maximum:", "Mean:      "
  )

  # construct abundance counts, same order as Fasta reads
  num_reads <- dataset$get_num_unique_seqs()
  total <- dataset$get_num_seqs()

  if (total == 0) {
    message <- "[ERROR]: dataset is empty. summary_seqs requires at least one
    sequence."
    cli::cli_abort(message)
  }

  fasta_results <- summarize_fasta(
    dataset$get_summary_report(),
    dataset$get_seqs_abunds(), processors
  )

  rownames(fasta_results) <- results_row_names

  results <- list(fasta_summary = fasta_results, contigs_summary = c())

  # if you have summary results to print
  if ((length(fasta_results$starts) == 8) && (!silent)) {
    cat("\t\tstart\tend\tlength\tambigs\tpolymer\tnum_ns\tnumseqs\n")
    cat(
      paste(results_row_names[1], fasta_results$starts[1],
        fasta_results$ends[1],
        fasta_results$nbases[1], fasta_results$ambigs[1],
        fasta_results$polymers[1], fasta_results$numns[1],
        fasta_results$numseqs[1],
        sep = "\t"
      ),
      "\n"
    )
    cat(
      paste(results_row_names[2], fasta_results$starts[2],
        fasta_results$ends[2],
        fasta_results$nbases[2], fasta_results$ambigs[2],
        fasta_results$polymers[2], fasta_results$numns[2],
        fasta_results$numseqs[2],
        sep = "\t"
      ),
      "\n"
    )
    cat(
      paste(results_row_names[3], fasta_results$starts[3],
        fasta_results$ends[3],
        fasta_results$nbases[3], fasta_results$ambigs[3],
        fasta_results$polymers[3], fasta_results$numns[3],
        fasta_results$numseqs[3],
        sep = "\t"
      ),
      "\n"
    )
    cat(
      paste(results_row_names[4], fasta_results$starts[4],
        fasta_results$ends[4],
        fasta_results$nbases[4], fasta_results$ambigs[4],
        fasta_results$polymers[4], fasta_results$numns[4],
        fasta_results$numseqs[4],
        sep = "\t"
      ),
      "\n"
    )
    cat(
      paste(results_row_names[5], fasta_results$starts[5],
        fasta_results$ends[5],
        fasta_results$nbases[5], fasta_results$ambigs[5],
        fasta_results$polymers[5], fasta_results$numns[5],
        fasta_results$numseqs[5],
        sep = "\t"
      ),
      "\n"
    )
    cat(
      paste(results_row_names[6], fasta_results$starts[6],
        fasta_results$ends[6],
        fasta_results$nbases[6], fasta_results$ambigs[6],
        fasta_results$polymers[6], fasta_results$numns[6],
        fasta_results$numseqs[6],
        sep = "\t"
      ),
      "\n"
    )
    cat(
      paste(results_row_names[7], fasta_results$starts[7],
        fasta_results$ends[7],
        fasta_results$nbases[7], fasta_results$ambigs[7],
        fasta_results$polymers[7], fasta_results$numns[7],
        fasta_results$numseqs[7],
        sep = "\t"
      ),
      "\n"
    )
    cat(paste(
      results_row_names[8], fasta_results$starts[8], fasta_results$ends[8],
      fasta_results$nbases[8], fasta_results$ambigs[8],
      fasta_results$polymers[8], fasta_results$numns[8],
      sep = "\t"
    ), "\n")
    cat("Unique seqs:\t", num_reads, "\n")
    cat("Total seqs:\t", total, "\n")
  }

  return(results)
}

# ============================================================================
