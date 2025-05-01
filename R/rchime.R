#' @title rchime
#'
#' @description
#' The rchime function allows you to remove chimeras from your data
#' using a denovo approach or alternatively a reference model.
#'
#' Our preferred way of doing this is to use the abundant sequences
#' as our reference (denovo). In addition, if a sequence is flagged as chimeric
#' in one sample, the default (dereplicate=FALSE) is to remove it from
#' all samples. Our experience suggests that this is a bit aggressive
#' since we’ve seen rare sequences get flagged as chimeric when they’re
#' the most abundant sequence in another sample. For a more conservative
#' approach, set (dereplicate=TRUE) which will only remove sequences
#' from the samples in which they are flagged as chimeric.
#'
#'  This function relies on a rewrite of Robert C. Edgar's
#'  \href{http://drive5.com/usearch/manual/uchime_algo.html}{uchime} source code
#'  which was donated to the public domain.
#'
#' @references
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
#' @param chimealns bool, Default = FALSE. When TRUE a DataFrame is returned
#'                  containing multiple alignments of query sequences to parents
#'                  in human-readable format. Alignments show columns with
#'                  differences that support or contradict a chimeric model.
#' @param rchime_options List, You can fine tune the uchime specific options
#' using  the [[rchime_options()]] function. Default = NULL.
#'
#' @seealso [uchimeout_help()] for a description of 'uchimeout'
#' @seealso [rchime_options()] to set uchime specific parameters.
#'
#' @return List containing accnos, uchimeout and optional chimealns.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @examples
#'
#' # Let's create a dataset from the example files
#'
#' dataset <- sequence_data_vector$new(
#'   filename = rchime_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = rchime_example("test.count_table")
#' )
#'
#' # detect and remove chimeras from the dataset using denovo approach
#' # (recommended)
#'
#' chimera_report <- rchime(dataset, dereplicate = TRUE)
#'
#' # or detect and remove chimeras using a reference
#'
#' dataset <- sequence_data_vector$new(
#'   filename = rchime_example("test.fasta")
#' )
#' dataset$set_group_assignments(
#'   filename = rchime_example("test.count_table")
#' )
#'
#' chimera_report <- rchime(dataset,
#'   reference = rchime_example("reference.fasta")
#' )
#'
#' @import cli
#' @importFrom parallelly, availableCores
#'
#' @export
rchime <- function(dataset = NULL, fasta = NULL, count = NULL,
                   reference = NULL, dereplicate = FALSE,
                   processors = parallelly::availableCores(),
                   silent = FALSE, chimealns = FALSE,
                   rchime_options = NULL) {
  created_dataset <- FALSE
  start_time <- Sys.time()

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
    dataset <- sequence_data_vector$new(filename = fasta)
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
      "[ERROR]: rchime requires your dataset to be ",
      "aligned."
    )
    cli::cli_abort(message)
  }

  parameters <- list(
    processors = processors, silent = silent,
    dereplicate = dereplicate
  )
  if (chimealns) {
    parameters <- c(parameters, chimealns = TRUE)
  }

  # if the user sets uchime options then add to parameters
  if (!is.null(rchime_options)) {
    if (!is.null(rchime_options$abskew)) {
      parameters <- c(parameters, abskew = rchime_options$abskew)
    }
    if (!is.null(rchime_options$minh)) {
      parameters <- c(parameters, minh = rchime_options$minh)
    }
    if (!is.null(rchime_options$mindiv)) {
      parameters <- c(parameters, mindiv = rchime_options$mindiv)
    }
    if (!is.null(rchime_options$xn)) {
      parameters <- c(parameters, xn = rchime_options$xn)
    }
    if (!is.null(rchime_options$dn)) {
      parameters <- c(parameters, dn = rchime_options$dn)
    }
    if (!is.null(rchime_options$xa)) {
      parameters <- c(parameters, xa = rchime_options$xa)
    }
    if (!is.null(rchime_options$chunks)) {
      parameters <- c(parameters, chunks = rchime_options$chunks)
    }
    if (!is.null(rchime_options$minchunk)) {
      parameters <- c(parameters, minchunk = rchime_options$minchunk)
    }
    if (!is.null(rchime_options$idsmoothwindow)) {
      parameters <- c(parameters,
        idsmoothwindow = rchime_options$idsmoothwindow
      )
    }
    if (!is.null(rchime_options$maxp)) {
      parameters <- c(parameters, maxp = rchime_options$maxp)
    }
    if (!is.null(rchime_options$skipgaps)) {
      parameters <- c(parameters, rchime_options$skipgaps)
    }
    if (!is.null(rchime_options$skipgaps2)) {
      parameters <- c(parameters, rchime_options$skipgaps2)
    }
    if (!is.null(rchime_options$minlen)) {
      parameters <- c(parameters, minlen = rchime_options$minlen)
    }
    if (!is.null(rchime_options$maxlen)) {
      parameters <- c(parameters, maxlen = rchime_options$maxlen)
    }
  }

  # build reference dataset if provided
  chimera_report <- NULL
  if (!is.null(reference)) {
    directory <- dirname(reference)
    filename <- basename(reference)

    reference_file <- file.path(file.path(directory), reference)
    df <- microseq::readFasta(reference_file)

    names <- unlist(lapply(
      df$Header,
      (function(x) {
        extract_name(x)
      })
    ))

    chimera_report <- rchimeReference(dataset, names, df$Sequence, parameters)
  } else {
    chimera_report <- rchimeDenovo(dataset, parameters)
  }

  if (created_dataset) {
    chimera_report[["dataset"]] <- dataset
  }

  if (!silent) {
    after_count <- dataset$get_num_seqs()
    if (after_count < num_seqs) {
      message <- paste(
        "rchime removed {.var {num_seqs-after_count}}",
        "chimeras from your dataset"
      )
    } else if (after_count == num_seqs) {
      message <- ("rchime complete, no chimeras found")
    }
    cli::cli_alert(message)
    message <- paste(
      "It took {.var {Sys.time() - start_time}}",
      "seconds to remove the chimeras"
    )
    cli::cli_alert(message)
  }

  return(chimera_report)
}
