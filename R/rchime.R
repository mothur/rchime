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
#' object containing reference sequence data.
#' @param dereplicate, Boolean, The dereplicate option allows you to remove
#'  chimeras by sample. For example, if dereplicate parameter is FALSE, then if
#'  one group finds the sequence to be chimeric, it will be removed from all
#'  groups. If dereplicate is set to TRUE, sequences found to be chimeric are
#'  only removed from the sample they are found to be chimeric in.
#'   Default = FALSE.
#' @param processors Integer, number of cores to use. Default = all available
#' @param remove_chimeras Boolean, remove chimeras from dataset. Default = TRUE.
#' @param silent Boolean, suppress concole outputs. Default = FALSE.
#' @param rchime_options List, You can fine tune the vsearch specific options
#' using  the [[rchime_options()]] function. Default = NULL.
#'
#' @seealso [rchime_options()] to set vsearch specific parameters.
#'
#' @return data.frame() containing chimera report.
#' The \href{https://mothur.org/strollur/}{strollur} dataset object will also be
#' updated. The sequences flagged as chimeric are removed and the chimera report
#' is added.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @examples
#'
#' # Let's create a dataset named rchime example
#'
#' data_denovo <- strollur::new_dataset("rchime denovo example")
#'
#' # Read a FASTA file and a mothur formatted count file
#'
#' fasta_data <- strollur::read_fasta(rchime_example("miseq.ng.fasta"))
#' abundance_data <- strollur::read_mothur_count(
#'   rchime_example("miseq.count_table")
#' )
#'
#' # Add the sequence and abundance data to the dataset
#'
#' strollur::add(data_denovo, table = fasta_data, type = "sequences")
#' strollur::assign(data_denovo,
#'   table = abundance_data,
#'   type = "sequence_abundance"
#' )
#'
#' # Detect and remove chimeras from the dataset using denovo approach
#' # (recommended)
#'
#' chimera_report <- rchime(data_denovo, dereplicate = TRUE)
#' data_denovo
#'
#' # Alternatively you can detect and remove chimeras using a reference
#'
#' data_reference <- strollur::new_dataset("rchime reference example")
#'
#' strollur::add(data_reference, table = fasta_data, type = "sequences")
#' strollur::assign(data_reference,
#'   table = abundance_data,
#'   type = "sequence_abundance"
#' )
#'
#' reference <- new_dataset("Silva V4 Region")
#'
#' reference_data <- strollur::read_fasta(rchime_example(
#'   "silva.v4.ng.fasta.gz"
#' ))
#'
#' strollur::add(reference, table = reference_data, type = "sequences")
#'
#' chimera_report <- rchime(data_reference, reference = reference)
#' data_reference
#'
#' @import cli
#' @importFrom parallelly, availableCores
#'
#' @export
rchime <- function(data, reference = NULL, dereplicate = FALSE,
                   processors = parallelly::availableCores(),
                   silent = FALSE, remove_chimeras = TRUE,
                   rchime_options = NULL) {
  if (!("dataset" %in% class(data))) {
    stop("data must be a stroller::dataset object")
  }

  start_time <- Sys.time()

  num_seqs <- strollur::count(data, type = "sequences")

  parameters <- list(
    processors = processors,
    dereplicate = dereplicate
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
  }

  num_samples <- strollur::count(data, type = "samples")
  # inputs needed ->
  ## reference -> names, seqs abunds, refnames, refseqs
  ## denovo no groups -> names, seqs, abunds

  ## denovo with groups -> names, seqs and abunds parsed by sample


  results <- NULL

  if (!is.null(reference)) {
    if (!("dataset" %in% class(reference))) {
      stop("reference must be a stroller::dataset object")
    }

    # reference -> names, seqs abunds, refnames, refseqs
    results <- rchimeReference(
      strollur::names(data, type = "sequences"),
      strollur::xdev_get_sequences(data),
      strollur::xdev_abundance(data)$abundances,
      strollur::names(reference, type = "sequences"),
      strollur::xdev_get_sequences(reference),
      parameters
    )
  } else {
    # chimera_report <- rchimeDenovo(dataset, parameters)

    if (num_samples == 0) {
      # denovo no groups -> names, seqs, abunds
      results <- rchimeDenovoSingleSample(
        strollur::names(data, type = "sequences"),
        strollur::xdev_get_sequences(data),
        strollur::xdev_abundance(data)$abundances,
        parameters
      )
    } else {
      # denovo with groups -> names, seqs and abunds parsed by sample
      results <- rchimeDenovo(
        strollur::xdev_get_by_sample(data),
        strollur::xdev_get_by_sample(
          data,
          type == "sequences"
        ),
        strollur::xdev_get_abundances_by_sample(data),
        parameters
      )
    }
  }

  # remove chimeras
  if (remove_chimeras) {
    # only remove chimeras from the samples they were flagged in
    if ((dereplicate) && (num_samples != 0)) {
      # set abundances parsed by sample
      strollur::xdev_set_abundances(
        data,
        names(results$set_abundance_values),
        unname(results$set_abundance_values),
        "chimeras_rchime"
      )
    } else {
      # remove all sequences flagged as chimeric
      if (length(results$accnos) > 0) {
        strollur::xdev_remove_sequences(
          data,
          results$accnos,
          rep(
            "chimeras_rchime",
            length(results$accnos)
          )
        )
      }
    }
  }

  if (!silent) {
    if (remove_chimeras) {
      after_count <- strollur::count(data, type = "sequences")
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
      cli::cli_alert(paste0("It took {.var {timing}} seconds to detect and ",
                            "remove the chimeras."))
    } else {
      timing <- difftime(Sys.time(), start_time, units = "secs")[[1]]
      cli::cli_alert("It took {.var {timing}} seconds to detect the chimeras.")
    }
  }

  results
}
