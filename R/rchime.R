#' @import RcppThread
#' @import cli
#'
#' @name rchime
#' @rdname rchime
#' @title Detect and remove chimeras from your
#'   \href{https://mothur.org/strollur/}{strollur} object or data.frame
#'   using a denovo approach or alternatively a reference model.
#' @description
#' The `rchime()` function allows you to detect and remove chimeras from your
#' data using a denovo approach or alternatively a reference model.
#'
#' Our preferred way of doing this is to use the abundant sequences as our
#' reference (denovo).
#'
#' This function uses code from the
#'  \href{https://github.com/torognes/vsearch}{vsearch} tools.
#'
#' @references
#'   Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
#'   versatile open source tool for metagenomics. PeerJ 4:e2584.
#'   doi: 10.7717/peerj.2584
#'
#' @references
#'  Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
#'  UCHIME improves sensitivity and speed of chimera detection.
#'  Bioinformatics 27:2194.
#'
#' @param data a \href{https://mothur.org/strollur/}{strollur} dataset object or
#'   a data.frame containing your sequence data.
#' @param reference a \href{https://mothur.org/strollur/}{strollur} dataset
#'   object or a data.frame containing reference sequence data.
#' @param dereplicate logical. The dereplicate option allows you to remove
#'   chimeras by sample. When `dereplicate=FALSE`, if a sequence is flagged as
#'   chimeric in one sample, it is removed from all samples. Our experience
#'   suggests that this is a bit aggressive since we’ve seen rare sequences get
#'   flagged as chimeric when they’re the most abundant sequence in another
#'   sample. For a more conservative approach, we recommend using the default
#'   `dereplicate=TRUE` which will only remove sequences from the samples in
#'   which they are flagged as chimeric.
#' @param remove_chimeras Boolean, remove chimeras from dataset. Default = TRUE.
#'   Only used when `data` is a strollur object.
#' @param silent Boolean, suppress console outputs. Default = FALSE.
#' @param rchime_options List, You can fine tune the vsearch specific options
#' using  the [[rchime_options()]] function. Default = NULL.
#' @param table_names, named list used to indicate the names of the columns in
#' the data.frame. Only used when `data` or `reference` are a data.frames.
#' @seealso [rchime_options()] to set vsearch specific parameters.
#'
#' @examples
#'
#' # Let's use a strollur object with 100 sequences
#'
#' data <- strollur::load_dataset(
#'   rchime_example("strollur_multi_sample_tiny.rds")
#' )
#'
#' chimera_report <- rchime(data)
#' data
#'
#' @return list() containing a chimera report, and vector of the chimeric
#' sequence's names.
#'
#' The \href{https://mothur.org/strollur/}{strollur} dataset object will also be
#' updated. The sequences flagged as chimeric are removed and the chimera report
#' is added.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @export
rchime <- function(data, reference = NULL, dereplicate = TRUE,
                   silent = FALSE, remove_chimeras = TRUE,
                   rchime_options = NULL,
                   table_names = list(
                     sequence_name = "sequence_name",
                     sequence = "sequence",
                     abundance = "abundance",
                     sample = "sample"
                   )) {
  UseMethod("rchime", data)
}


#' @title Detect and remove chimeras from your
#' \href{https://mothur.org/strollur/}{strollur} dataset object.
#' @name rchime.strollur
#' @rdname rchime.strollur
#' @description
#' The `rchime()` function allows you to detect and remove chimeras from your
#' data using a denovo approach or alternatively a reference model.
#'
#' Our preferred way of doing this is to use the abundant sequences as our
#' reference (denovo).
#'
#' This function uses code from the
#'  \href{https://github.com/torognes/vsearch}{vsearch} tools.
#'
#' @references
#'   Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
#'   versatile open source tool for metagenomics. PeerJ 4:e2584.
#'   doi: 10.7717/peerj.2584
#'
#' @references
#'  Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
#'  UCHIME improves sensitivity and speed of chimera detection.
#'  Bioinformatics 27:2194.
#'
#' @param data a \href{https://mothur.org/strollur/}{strollur} dataset object
#' containing your sequence data.
#' @param reference a \href{https://mothur.org/strollur/}{strollur} dataset
#' object or a data.frame containing reference sequence data.
#' @param dereplicate logical. The dereplicate option allows you to remove
#'  chimeras by sample. When `dereplicate=FALSE`, if a sequence is flagged as
#'   chimeric in one sample, it is removed from all samples. Our experience
#'   suggests that this is a bit aggressive since we’ve seen rare sequences get
#'   flagged as chimeric when they’re the most abundant sequence in another
#'   sample. For a more conservative approach, we recommend using the default
#'   `dereplicate=TRUE` which will only remove sequences from the samples in
#'   which they are flagged as chimeric.
#' @param remove_chimeras Boolean, remove chimeras from dataset. Default = TRUE.
#' @param silent Boolean, suppress console outputs. Default = FALSE.
#' @param rchime_options List, You can fine tune the vsearch specific options
#' using  the [[rchime_options()]] function. Default = NULL.
#' @seealso [rchime_options()] to set vsearch specific parameters.
#' @param table_names, Only used when `reference` is a data.frame. By default:
#'
#' table_names <- list(sequence_name = "sequence_name",
#'                     sequence = "sequence")
#'
#' In table_names, 'sequence_name' is a string containing the name of the column
#' in 'table' that contains the sequence names. Default column name is
#' 'sequence_name'.
#'
#' In table_names, 'sequence' is a string containing the name of the
#' column in 'table' that contains the sequences. Default column name is
#'  'sequence'.
#'
#' @return list() containing a chimera report, and vector of the chimeric
#' sequence's names.
#'
#' The \href{https://mothur.org/strollur/}{strollur} dataset object will also be
#' updated. The sequences flagged as chimeric are removed and the chimera report
#' is added.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @examples
#'
#' # Let's load a strollur object with 500 sequences
#'
#' data <- strollur::load_dataset(
#'   rchime_example("strollur_multi_sample_small.rds")
#' )
#'
#' # Detect and remove chimeras from the dataset using denovo approach by sample
#' # (recommended)
#'
#' chimera_report <- rchime(data)
#' data
#'
#' @import cli
#' @export
rchime.strollur <- function(data, reference = NULL, dereplicate = TRUE,
                            silent = FALSE, remove_chimeras = TRUE,
                            rchime_options = NULL,
                            table_names = list(
                              sequence_name = "sequence_name",
                              sequence = "sequence"
                            )) {
  if (!("strollur" %in% class(data))) {
    stop("data must be a stroller::strollur object")
  }

  start_time <- Sys.time()

  num_seqs <- strollur::count(data, type = "sequence")

  parameters <- list(
    dereplicate = dereplicate,
    processors = parallelly::availableCores()
  )

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
    if (!is.null(rchime_options$maxp)) {
      parameters <- c(parameters, maxp = rchime_options$maxp)
    }
    if (!is.null(rchime_options$processors)) {
      parameters[["processors"]] <- rchime_options$processors
      processors <- rchime_options$processors
    }
    if (!is.null(rchime_options$dereplicate)) {
      parameters[["dereplicate"]] <- rchime_options$dereplicate
      dereplicate <- rchime_options$dereplicate
    }
  }

  num_samples <- strollur::count(data, type = "sample")

  results <- NULL

  if (!is.null(reference)) {
    if (!("strollur" %in% class(reference)) &&
          !("data.frame" %in% class(reference))) {
      stop("reference must be a strollur object or data.frame")
    }

    dereplicate <- FALSE

    # reference -> names, seqs abunds, refnames, refseqs
    if ("strollur" %in% class(reference)) {
      results <- rchimeReference(
        strollur::names(data, type = "sequence"),
        strollur::xdev_get_sequences(data, degap = TRUE),
        strollur::names(reference, type = "sequence"),
        strollur::xdev_get_sequences(reference, degap = TRUE),
        parameters
      )
    } else {
      default_tn <- list(
        sequence_name = "sequence_name",
        sequence = "sequence"
      )

      table_names <- modifyList(default_tn, table_names)

      results <- rchimeReference(
        strollur::names(data, type = "sequence"),
        strollur::xdev_get_sequences(data, degap = TRUE),
        fill_required_parameters(
          reference,
          table_names[["sequence_name"]]
        ),
        gsub("[.-]", "", fill_required_parameters(
          reference,
          table_names[["sequence"]]
        )),
        parameters
      )
    }
  } else {
    if (!silent) {
      message <- "The denovo method runs with a single processor.\n\n"
      cli::cli_alert_info(message)
    }
    if (num_samples == 0) {
      # denovo no groups -> names, seqs, abunds
      results <- rchimeDenovoSingleSample(
        strollur::names(data, type = "sequence"),
        strollur::xdev_get_sequences(data, degap = TRUE),
        strollur::xdev_abundance(data)$abundance,
        parameters
      )
    } else {
      # denovo with groups -> names, seqs and abunds parsed by sample
      results <- rchimeDenovo(
        strollur::xdev_get_by_sample(data),
        strollur::xdev_get_by_sample(
          data,
          type = "sequence",
          degap = TRUE
        ),
        strollur::xdev_get_abundances_by_sample(data),
        parameters
      )
    }
  }

  # add chimera report to data
  strollur::add(data,
    table = results$chimera_report,
    type = "report",
    report_type = "chimera_report",
    table_names = list(sequence_name = "Query"),
    verbose = !silent
  )

  # remove chimeras
  if (remove_chimeras) {
    # only remove chimeras from the samples they were flagged in
    if ((dereplicate) && (num_samples != 0)) {
      # set abundances parsed by sample
      strollur::xdev_set_abundances(
        data,
        results$set_abundance_values$sequence_name,
        results$set_abundance_values$abundance,
        "rchime-chimeras"
      )
      results[["set_abundance_values"]] <- NULL
    } else {
      # remove all sequences flagged as chimeric
      if (length(results$chimeras) > 0) {
        strollur::xdev_remove_sequences(
          data,
          results$chimeras,
          rep(
            "chimeras_rchime",
            length(results$chimeras)
          )
        )
      }
    }
  }

  if (!silent) {
    if (remove_chimeras) {
      after_count <- strollur::count(data, type = "sequence")
      if (after_count < num_seqs) {
        message <- paste(
          "rchime removed {.var {num_seqs-after_count}}",
          "chimeras from your dataset."
        )
      } else if (after_count == num_seqs) {
        message <- ("rchime complete, no chimeras found.")
      }
      cli::cli_alert(message)
      timing <- difftime(Sys.time(), start_time, units = "secs")[[1]]
      cli::cli_alert(paste0(
        "It took {.var {timing}} seconds to detect and ",
        "remove the chimeras."
      ))
    } else {
      # report detected
      if ((dereplicate) && (num_samples != 0)) {
        after_count <- sum(unlist(results$set_abundance_values$abundances))
      } else {
        # number of chimeras removed
        df <- strollur::abundance(data, type = "sequence")
        num_chimeras <- sum(df$abundance[
          df$sequence_name %in% results$chimeras
        ])

        after_count <- num_seqs - num_chimeras
      }
      if (after_count < num_seqs) {
        message <- paste(
          "rchime detected {.var {num_seqs-after_count}}",
          "chimeras in your dataset."
        )
      } else if (after_count == num_seqs) {
        message <- ("rchime complete, no chimeras found.")
      }
      cli::cli_alert(message)
      timing <- difftime(Sys.time(), start_time, units = "secs")[[1]]
      cli::cli_alert("It took {.var {timing}} seconds to detect the chimeras.")
    }
  }

  results
}

# =============================================================================#
#' @title Detect chimeras in your data.frames.
#' @name rchime.data.frame
#' @rdname rchime.data.frame
#' @description
#' The `rchime()` function allows you to detect chimeras from your
#' data using a denovo approach or alternatively a reference model.
#'
#' Our preferred way of doing this is to use the abundant sequences as our
#' reference (denovo).
#'
#' This function uses code from the
#'  \href{https://github.com/torognes/vsearch}{vsearch} tools.
#'
#' @references
#'   Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
#'   versatile open source tool for metagenomics. PeerJ 4:e2584.
#'   doi: 10.7717/peerj.2584
#'
#' @references
#'  Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
#'  UCHIME improves sensitivity and speed of chimera detection.
#'  Bioinformatics 27:2194.
#'
#' @param data a data.frame containing your sequence data.
#' @param reference a data.frame containing reference sequence data.
#' @param dereplicate logical. The dereplicate option allows you to flag
#'  chimeras by sample. When `dereplicate=FALSE`, if a sequence is flagged as
#'   chimeric in one sample, it is removed from all samples. Our experience
#'   suggests that this is a bit aggressive since we’ve seen rare sequences get
#'   flagged as chimeric when they’re the most abundant sequence in another
#'   sample. For a more conservative approach, we recommend using the default
#'   `dereplicate=TRUE` which will only remove sequences from the samples in
#'   which they are flagged as chimeric.
#' @param silent Boolean, suppress console outputs. Default = FALSE.
#' @param remove_chimeras Only used when `data` is a strollur object.
#' @param rchime_options List, You can fine tune the vsearch specific options
#' using  the [[rchime_options()]] function. Default = NULL.
#' @param table_names, named list used to indicate the names of the columns in
#' the data.frame. Only used when `data` is a data.frame. By default:
#'
#' table_names <- list(sequence_name = "sequence_name",
#'                     sequence = "sequence"
#'                     abundance = "abundance",
#'                     sample = "sample")
#'
#' In table_names, 'sequence_name' is a string containing the name of the column
#' in 'table' that contains the sequence names. Default column name is
#' 'sequence_name'.
#'
#' In table_names, 'sequence' is a string containing the name of the
#' column in 'table' that contains the sequences. Default column name is
#'  'sequence'.
#'
#' In table_names, 'abundance' is a string containing the name of the column in
#' 'table' that contains the abundances. Default column name is 'abundance'.
#'
#' In table_names, 'sample' is a string containing the name of the column in
#' 'table' that contains the samples. Default column name is 'sample'.
#'
#' @seealso [rchime_options()] to set vsearch specific parameters.
#'
#' @return list() containing a chimera report, and vector of the chimeric
#'   sequence's names. If running with multiple samples and dereplicate = TRUE,
#'   then a table containing the modified sequence abundances will also be
#'   returned.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @examples
#'
#' # Detect chimeras from the dataset using denovo approach by sample
#' # (recommended)
#'
#' data <- readRDS(rchime_example("miseq_data_frame_by_sample_small.rds"))
#'
#' chimera_report <- rchime(data)
#'
#' @import cli
#' @export
rchime.data.frame <- function(data, reference = NULL, dereplicate = TRUE,
                              silent = FALSE, remove_chimeras = NULL,
                              rchime_options = NULL,
                              table_names = list(
                                sequence_name = "sequence_name",
                                sequence = "sequence",
                                abundance = "abundance",
                                sample = "sample"
                              )) {
  if (!("data.frame" %in% class(data))) {
    stop("data must be a data.frame object")
  }

  start_time <- Sys.time()

  default_tn <- list(
    sequence_name = "sequence_name",
    abundance = "abundance",
    sample = "sample",
    sequence = "sequence"
  )

  table_names <- modifyList(default_tn, table_names)

  parameters <- list(
    dereplicate = dereplicate,
    processors = parallelly::availableCores()
  )

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
    if (!is.null(rchime_options$maxp)) {
      parameters <- c(parameters, maxp = rchime_options$maxp)
    }
    if (!is.null(rchime_options$processors)) {
      parameters[["processors"]] <- rchime_options$processors
      processors <- rchime_options$processors
    }
    if (!is.null(rchime_options$dereplicate)) {
      parameters[["dereplicate"]] <- rchime_options$dereplicate
      dereplicate <- rchime_options$dereplicate
    }
  }

  num_samples <- 0

  # check for sample info
  if (table_names[["sample"]] %in% names(data)) {
    # create strollur object to read and parse inputs
    strollur_data <- strollur::new_dataset()

    # add sequence data
    strollur::add(strollur_data,
      table = unique(data[, c(
        table_names[["sequence_name"]],
        table_names[["sequence"]]
      )]),
      type = "sequence",
      table_names = table_names, verbose = FALSE
    )

    # assign sequence abundance
    strollur::assign(strollur_data,
      table = data,
      type = "sequence_abundance",
      table_names = table_names, verbose = FALSE
    )

    num_samples <- strollur::count(strollur_data, type = "sample")
  }

  results <- NULL

  if (!is.null(reference)) {
    if (!("data.frame" %in% class(data))) {
      stop("reference must be a data.frame")
    }

    dereplicate <- FALSE

    # reference -> names, seqs, refnames, refseqs
    results <- rchimeReference(
      fill_required_parameters(
        data,
        table_names[["sequence_name"]]
      ),
      gsub("[.-]", "", fill_required_parameters(
        data,
        table_names[["sequence"]]
      )),
      fill_required_parameters(
        reference,
        table_names[["sequence_name"]]
      ),
      gsub("[.-]", "", fill_required_parameters(
        reference,
        table_names[["sequence"]]
      )),
      parameters
    )
  } else {
    if (!silent) {
      message <- "The denovo method runs with a single processor.\n\n"
      cli::cli_alert_info(message)
    }
    if (num_samples == 0) {
      # denovo no groups -> names, seqs, abunds
      results <- rchimeDenovoSingleSample(
        fill_required_parameters(data, table_names[["sequence_name"]]),
        gsub("[.-]", "", fill_required_parameters(
          data,
          table_names[["sequence"]]
        )),
        fill_required_parameters(data, table_names[["abundance"]]),
        parameters
      )
    } else {
      # denovo with groups -> names, seqs and abunds parsed by sample
      results <- rchimeDenovo(
        strollur::xdev_get_by_sample(strollur_data),
        strollur::xdev_get_by_sample(
          strollur_data,
          type = "sequence",
          degap = TRUE
        ),
        strollur::xdev_get_abundances_by_sample(strollur_data),
        parameters
      )

      if (dereplicate) {
        # add sample names to updated abundance info
        sample_names <- strollur::names(strollur_data, type = "sample")
        results$set_abundance_values[["samples"]] <- sample_names
      }
    }
  }

  if (!silent) {
    # report detected
    if ((dereplicate) && (num_samples != 0)) {
      after_count <- sum(unlist(results$set_abundance_values$abundance))
      num_seqs <- strollur::count(strollur_data, type = "sequence")
    } else {
      # number of chimeras removed
      if (table_names[["abundance"]] %in% names(data)) {
        num_chimeras <- sum(data[[table_names[["abundance"]]]][
          data[[table_names[["sequence_name"]]]] %in% results$chimeras
        ])
        num_seqs <- sum(data[[table_names[["abundance"]]]])
      } else {
        # assume unique
        num_chimeras <- length(results$chimeras)
        num_seqs <- nrow(results$chimera_report)
      }

      after_count <- num_seqs - num_chimeras
    }
    if (after_count < num_seqs) {
      message <- paste(
        "rchime detected {.var {num_seqs-after_count}}",
        "chimeras in your dataset."
      )
    } else if (after_count == num_seqs) {
      message <- ("rchime complete, no chimeras found.")
    }
    cli::cli_alert(message)
    timing <- difftime(Sys.time(), start_time, units = "secs")[[1]]
    cli::cli_alert("It took {.var {timing}} seconds to detect the chimeras.")
  }

  results
}
# =============================================================================#
#' @export
rchime.default <- function(data, reference = NULL, dereplicate = TRUE,
                           silent = FALSE, remove_chimeras = TRUE,
                           rchime_options = NULL, table_names = list(
                             sequence_name = "sequence_name",
                             sequence = "sequence",
                             abundance = "abundance",
                             sample = "sample"
                           )) {
  stop("data must be a strollur object or a data.frame")
}
# =============================================================================#
