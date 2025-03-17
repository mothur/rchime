#' @title chimera_uchime
#'
#' @description
#' The chimera_uchime function detects and removes chimeras from your data.
#'
#'  This function uses modified
#'  \href{http://drive5.com/usearch/manual/uchime_algo.html}{uchime} source code
#'   written by Robert C. Edgar and donated to the public domain.
#'
#'  Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
#'  UCHIME improves sensitivity and speed of chimera detection.
#'  Bioinformatics 27:2194.
#'
#' @param dataset R6 'sequence_dataset' object containing sequence data
#' @param fasta String, name of file containing FASTA reads
#' @param count String, name of mothur formatted
#'    \href{https://mothur.org/wiki/count_file/}{count} file.
#' @param reference String, name of optional FASTA reference file.
#' @param dereplicate, Boolean, The dereplicate option allows you to remove
#'  chimeras by sample. For example, if dereplicate parameter is false, then if
#'  one group finds the sequence to be chimeric, it will be removed from all
#'  groups. If dereplicate is set to true, sequences found to be chimeric are
#'  only removed from the sample they are found to be chimeric in.
#'   Default = FALSE.
#' @param processors Integer, number of cores to use. Default = all available
#' @param silent Boolean, suppress progress updates and console output
#'
#' @return DataFrame containing chimera results and an updated dataset with the
#'  chimeras removed
#'
#'  Column descriptions from uchime's manual:
#'
#'  1	Score	Higher score means more strongly chimeric alignment.
#'  2	Q	Query Name.
#'  3	A	Parent A Name.
#'  4	B	Parent B Name.
#'  5	T	Top parent (T) Name. This is the closest reference sequence; usually
#'          either A or B.
#'  6	IdQM	Percent identity of query and the model (M) constructed as a
#'              segment of A and a segment of B.
#'  7	IdQA	Percent identity of Q and A.
#'  8	IdQA	Percent identity of Q and B.
#'  9	IdAB	Percent identity of A and B
#'  10	IdQT	Percent identity of Q and T.
#'  11	LY	Yes votes in left segment.
#'  12	LN	No votes in left segment.
#'  13	LA	Abstain votes in left segment.
#'  14	RY	Yes votes in right segment.
#'  15	RN	No votes in right segment.
#'  16	RA	Abstain votes in right segment.
#'  17	Div	Divergence, defined as (IdQM - IdQT).
#'  18	YN	Y, N or ?, indicating whether the query was classified as
#'          chimeric (Y), not chimeric (N) or borderline case (?).
#'          The query is classified as chimeric if h >= threshold specified by
#'          -minh, Div > minimum divergence specified by ‑mindiv and the number
#'          of diffs ( (Y+N+A) in each segment (L and R) is greater than the
#'          minimum specified by -mindiffs. A query is unclassified if the
#'          maxh > h > minh, i.e. maxh is the maximum score for a non-chimera,
#'          and minh is the minimum score for a chimera; in between is
#'          unclassified.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @examples
#'
#' # Let's create a dataset from the example files
#'
#' dataset <- sequence_data_table$new(
#'   filename = mothur2_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = mothur2_example("test.count_table")
#' )
#'
#' # detect and remove chimeras from the dataset using denovo approach
#' # (recommended)
#'
#' chimera_report <- chimera_uchime(dataset, dereplicate = TRUE)
#'
#' # or detect and remove chimeras using a reference
#'
#' dataset <- sequence_data_table$new(
#'   filename = mothur2_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = mothur2_example("test.count_table")
#' )
#' reference <- sequence_data_table$new(
#'   filename = mothur2_example("reference.fasta")
#' )
#'
#' chimera_report <- chimera_uchime(dataset, reference = reference)
#'
#' @import cli
#' @importFrom parallelly, availableCores
#'
#' @export
chimera_uchime <- function(dataset = NULL, fasta = NULL, count = NULL,
                           reference = NULL, dereplicate = FALSE,
                           processors = parallelly::availableCores(),
                           silent = FALSE) {
  # no inputs provided
  if (is.null(dataset) && is.null(fasta)) {
    message <- paste(
      "[ERROR]: you must provide a 'sequence_dataset'",
      "object or a fasta file."
    )
    cli::cli_abort(message)
    # provided dataset and fasta file, which to use
  } else if (!is.null(dataset) && !is.null(fasta)) {
    message <- paste("[ERROR]: you can use the 'dataset' or 'fasta' option",
      ", not both.",
      sep = ""
    )
    cli::cli_abort(message)
    # create dataset from fasta file and count file if provided
  } else if (is.null(dataset) && !is.null(fasta)) {
    dataset <- sequence_data_table$new(filename = fasta)
    if (!is.null(count)) {
      dataset$set_group_assignments(filename = count)
    }
    if (dataset$get_num_seqs() == 0) {
      message <- paste(
        "[ERROR]: empty dataset created from",
        "{.var {fasta}}, please correct."
      )
      cli::cli_abort(message)
    }
    # check dataset is not empty
  } else {
    if (dataset$get_num_seqs() == 0) {
      message <- paste("[ERROR]: your dataset is empty, please correct.")
      cli::cli_abort(message)
    }
  }

  num_seqs <- dataset$get_num_seqs()

  parameters <- list(
    processors = processors, silent = silent,
    dereplicate = dereplicate
  )

  # build reference dataset if provided
  chimera_report <- NULL
  if (!is.null(reference)) {
    reference <- sequence_data_table$new(filename = reference)
    chimera_report <- chimeraUchimeReference(dataset, reference, parameters)
  } else {
    chimera_report <- chimeraUchime(dataset, parameters)
  }

  if (!silent) {
    after_count <- dataset$get_num_seqs()
    if (after_count < num_seqs) {
      message <- paste(
        "chimera_uchime removed {.var {num_seqs-after_count}}",
        "chimeras from your dataset"
      )
    } else if (after_count == num_seqs) {
      message <- ("chimera_uchime complete, no chimeras found")
    }
    cli::cli_alert(message)
  }

  return(chimera_report)
}
