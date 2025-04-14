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
#'
#' ## Parameters used directly by uchime. It is recommended to use the default
#' values for the parameters below unless you have a pressing reason to change
#' them.
#'
#' @param chimealns bool, Default = FALSE. When TRUE a DataFrame is returned
#'                  containing multiple alignments of query sequences to parents
#'                  in human-readable format. Alignments show columns with
#'                  differences that support or contradict a chimeric model.
#' @param abskew Float, the minimum abundance skew (denovo only). Default = 2.0.
#'               abskew <- min (abund(parent1), abund(parent2)) / abund(query)
#' @param minh Float, Mininum score to report chimera. Default 0.3. Values from
#'     0.1 to 5 might be reasonable. Lower values increase sensitivity
#'    but may report more false positives. If you decrease --xn,
#'    you may need to increase --minh, and vice versa.
#' @param mindiv Float, Minimum divergence ratio, default 0.5. Div ratio is
#'    100%% %%identity between query sequence and the closest candidate for
#'    being a parent. If you don't care about very close chimeras,
#'    then you could increase --mindiv to, say, 1.0 or 2.0, and
#'    also decrease --min h, say to 0.1, to increase sensitivity.
#'    How well this works will depend on your data. Best is to
#'    tune parameters on a good benchmark.
#' @param xn Float, Weight of a no vote, also called the beta parameter.
#'    Default 8.0. Decreasing this weight to around 3 or 4 may give better
#'    performance on denoised data.
#' @param dn Float, Pseudo-count prior on number of no votes. Default 1.4.
#'    Probably no good reason to change this unless you can retune to a good
#'    benchmark for your data. Reasonable values are probably in the
#'    range from 0.2 to 2.
#' @param xa Float, Weight of an abstain vote. Default 1. So far, results do not
#'    seem to be very sensitive to this parameter, but if you have
#'    a good training set might be worth trying. Reasonable values
#'    might range from 0.1 to 2.
#' @param chunks Integer, Number of chunks to extract from the query sequence
#'    when searching for parents. Default 4.
#' @param minchunk Integer, Minimum length of a chunk. Default 64.
#' @param idsmoothwindow Integer, Length of id smoothing window. Default 32.
#' @param maxp Integer, Maximum number of candidate parents to consider.
#'   Default 2. In tests so far, increasing maxp gives only a very small
#'   improvement in sensivity but tends to increase the error rate quite a bit.
#' @param skipgaps, bool, If skipgaps is specified, columns containing gaps do
#'  not count as diffs. Default = TRUE.
#' @param skipgaps2 bool, If skipgaps2 is specified, if column is immediately
#'  adjacent to a column containing a gap, it is not counted as a diff.
#'  Default = TRUE.
#' @param minlen Integer, Minimum sequence length. Default = 10.
#' @param maxlen Integer, Maximum sequence length. Default = 10000.
#'
#' @seealso [uchimeout_help()] for a description of the DataFrame returned
#'
#' @return DataFrame containing chimera results and an updated dataset with the
#'  chimeras removed
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @examples
#'
#' # Let's create a dataset from the example files
#'
#' dataset <- sequence_data_table$new(
#'   filename = rchime_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = rchime_example("test.count_table")
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
#'   filename = rchime_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = rchime_example("test.count_table")
#' )
#'
#' chimera_report <- chimera_uchime(dataset,
#'   reference = rchime_example("reference.fasta")
#' )
#'
#' @import cli
#' @importFrom parallelly, availableCores
#'
#' @export
chimera_uchime <- function(dataset = NULL, fasta = NULL, count = NULL,
                           reference = NULL, dereplicate = FALSE,
                           processors = parallelly::availableCores(),
                           silent = FALSE, chimealns = FALSE, abskew = 2.0,
                           minh = 0.3, mindiv = 0.5, xn = 8.0,
                           dn = 1.4, xa = 1.0, chunks = 4, minchunk = 64,
                           idsmoothwindow = 32, maxp = 2, skipgaps = TRUE,
                           skipgaps2 = TRUE, minlen = 10, maxlen = 10000) {
  created_dataset <- FALSE

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

    created_dataset <- TRUE
  } else {
    if (dataset$get_num_seqs() == 0) {
      message <- paste("[ERROR]: your dataset is empty, please correct.")
      cli::cli_abort(message)
    }
  }

  num_seqs <- dataset$get_num_seqs()

  # make sure the dataset is aligned
  if (!dataset$is_aligned()$aligned) {
    message <- paste0(
      "[ERROR]: chimera_uchime requires your dataset to be ",
      "aligned."
    )
    cli::cli_abort(message)
  }

  parameters <- list(
    processors = processors, silent = silent,
    dereplicate = dereplicate
  )

  # if the user sets uchime options then add to parameters
  if (abskew != 1.9) {
    parameters <- c(parameters, abskew = abskew)
  }
  if (chimealns) {
    parameters <- c(parameters, chimealns = TRUE)
  }
  if (minh != 0.3) {
    parameters <- c(parameters, minh = minh)
  }
  if (mindiv != 0.5) {
    parameters <- c(parameters, mindiv = mindiv)
  }
  if (xn != 8.0) {
    parameters <- c(parameters, xn = xn)
  }
  if (dn != 1.4) {
    parameters <- c(parameters, dn = dn)
  }
  if (xa != 1.0) {
    parameters <- c(parameters, xa = xa)
  }
  if (chunks != 4) {
    parameters <- c(parameters, chunks = chunks)
  }
  if (minchunk != 64) {
    parameters <- c(parameters, minchunk = minchunk)
  }
  if (idsmoothwindow != 32) {
    parameters <- c(parameters, idsmoothwindow = idsmoothwindow)
  }
  if (maxp != 2) {
    parameters <- c(parameters, maxp = maxp)
  }
  if (!skipgaps) {
    parameters <- c(parameters, skipgaps = FALSE)
  }
  if (!skipgaps2) {
    parameters <- c(parameters, skipgaps2 = FALSE)
  }
  if (minlen != 10) {
    parameters <- c(parameters, minlen = minlen)
  }
  if (maxlen != 10000) {
    parameters <- c(parameters, maxlen = maxlen)
  }

  # build reference dataset if provided
  chimera_report <- NULL
  if (!is.null(reference)) {
    reference <- sequence_data_table$new(filename = reference)
    chimera_report <- chimeraUchimeReference(dataset, reference, parameters)
  } else {
    chimera_report <- chimeraUchime(dataset, parameters)
  }

  if (created_dataset) {
    chimera_report[["dataset"]] <- dataset
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
