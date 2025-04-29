#' @title sequence_dataset
#' @description 'sequence_dataset' is an R6 class that represents sequence FASTA
#'  and abundance data.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @import cli
#' @import stringi
#' @export
sequence_dataset <- R6Class("sequence_dataset",
  public = list(

    #' @description
    #' Create a new FASTA sequence dataset parent object
    #' @examples
    #'   dataset <- sequence_data_vector$new()
    #' @return A new `sequence_data_vector` object which inherits
    #'  from sequence_dataset
    initialize = function() {},

    #' @description
    #' Report name missing and abort
    #' @param name String, containing name of missing sequence
    abort_name_missing = function(name) {
      cli::cli_abort("[ERROR]: {.var {name}} is not in your dataset,
                            object mismatch.")
    },

    #' @description
    #' Report length mismatch and abort
    #' @param name1 String, containing name of first variable
    #' @param name2 String, containing name of second variable
    #' @param length1 Numeric, containing length of first variable
    #' @param length2 Numeric, containing length of second variable
    abort_length_mismatch = function(name1, name2, length1, length2) {
      cli::cli_abort("[ERROR]: The length of {.var {name1}}
          must equal the length of {.var {name2}}. {.var {name1}} has length
          {.var {length1}} but {.var {name2}} has length {.var {length2}}")
    },

    #' @description
    #' Report missing name
    #' @param name String, containing name of missing sequence
    #' @param return_value String, containing return_value
    alert_missing_name = function(name, return_value) {
      cli::cli_alert("{.var {name}} is not in your dataset,
                       returning {.var {return_value}}.")
    },

    #' @description
    #' Report eliminated name
    #' @param name String, containing name of eliminated sequence
    #' @param return_value String, containing return_value
    alert_eliminated_name = function(name, return_value) {
      cli::cli_alert("{.var {name}} has been eliminated from your
          dataset, returning {.var {return_value}}.")
    },

    #' @description
    #' Get number of ambiguous bases
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_ambigs = function(seq) {
      length <- stri_length(seq)
      num_non_ambig <- stringi::stri_count(seq, regex = "[*.*-ATGC]")
      return(length - num_non_ambig)
    },

    #' @description
    #' Get length of longest homopolymer
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_longest_polymer = function(seq) {
      longest_homop <- 1
      size <- 2

      # create homop vector of length 2 - start
      homops <- private$get_homop(size)

      # while you can still find homopolymers of this size
      while (any(stri_detect_fixed(seq, homops))) {
        longest_homop <- size
        size <- size + 1

        # increase homop length
        homops <- private$get_homop(size)
      }
      return(longest_homop)
    },

    #' @description
    #' Get number of bases in sequence
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_numbases = function(seq) {
      copy <- seq
      return(stringi::stri_length(
        stringi::stri_replace_all_regex(copy, "[.-]", "")
      ))
    },

    #' @description
    #' Get number of 'N' bases
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_numns = function(seq) {
      numns <- gregexpr("[N]", seq)[[1]][1]
      if (numns == -1) {
        return(0)
      }
      return(numns)
    },

    #' @description
    #' Get location of first base in sequence
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_start = function(seq) {
      return(regexpr("[A-Z]", seq)[1])
    },

    #' @description
    #' Get location of last base in sequence
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_end = function(seq) {
      return(stringi::stri_locate_last_regex(seq, "[A-Z]")[1])
    },

    #' @description
    #' Read a FASTA formatted sequence file
    #' @param filename FASTA file name (required)
    #' @param path path to FASTA file (optional)
    #' @examples
    #'   dataset <- sequence_dataset$new()
    #'   fasta_data <- dataset$read_fasta_file(filename =
    #'   rchime_example("test.fasta"))
    #' @return A data.frame containing the FASTA sequence data and summary
    read_fasta_file = function(filename, path = NULL) {
      directory <- dirname(filename)
      filename <- basename(filename)

      # if path not set, make directory of filename
      if (is.null(path)) {
        path <- directory
      }

      fastafile <- file.path(file.path(path), filename)
      df <- microseq::readFasta(fastafile)
      rows <- nrow(df)

      # extract name from comments
      names <- unlist(lapply(
        df$Header,
        (function(x) {
          private$extract_name(x)
        })
      ))

      # calc numNs
      numns <- unlist(lapply(
        df$Sequence,
        (function(x) {
          self$calc_numns(x)
        })
      ))

      # calc numBases
      lengths <- unlist(lapply(
        df$Sequence,
        (function(x) {
          self$calc_numbases(x)
        })
      ))

      # calc starts
      starts <- unlist(lapply(
        df$Sequence,
        (function(x) {
          self$calc_start(x)
        })
      ))

      # calc ends
      ends <- unlist(lapply(
        df$Sequence,
        (function(x) {
          self$calc_end(x)
        })
      ))


      ambigs <- unlist(lapply(
        df$Sequence,
        (function(x) self$calc_ambigs(x))
      ))

      polymers <- unlist(lapply(
        df$Sequence,
        (function(x) {
          self$calc_longest_polymer(x)
        })
      ))

      comments <- rep("", rows)
      trash <- rep("", rows)

      data <- data.frame(
        names = names, sequences = df$Sequence,
        comments = comments, trash_codes = trash,
        starts = starts, ends = ends, lengths = lengths,
        ambigs = ambigs, polymers = polymers,
        numns = numns
      )

      return(data)
    },

    #' @description
    #' Write dataset to FASTA file
    #' @param names vector of strings containing sequence names
    #' @param sequences vector of strings containing sequence nucleotides
    #' @param filename String, name of file you want to write FASTA
    #' sequences to
    #' @examples
    #'   dataset <- sequence_dataset$new()
    #'   fasta_data <- dataset$read_fasta_file(filename =
    #'   rchime_example("test.fasta"))
    #'  dataset$write_fasta_file(fasta_data$names,
    #'   fasta_data$sequences, filename = "output.fasta")
    #'  remove_file("output.fasta")
    write_fasta_file = function(names, sequences, filename) {
      if (length(names) != length(sequences)) {
        cli::cli_abort("[ERROR]: The length of names must equal the
                length of sequences. names has {.var {length(names)}} sequences
                but the number of sequences is {.var {length(sequences)}}")
      }

      if (file.exists(filename)) {
        file.remove(filename)
      }

      file_conn <- file(filename, open = "a")

      for (i in seq_along(names)) {
        print_name <- paste(">", names[i], sep = "")
        output <- paste(paste(print_name),
          sequences[i],
          sep = "\n"
        )

        writeLines(output, file_conn)
      }
      close(file_conn)
    }
  ),
  private = list(

    # extract name from line
    extract_name = function(name) {
      # name. comment
      # mySeq This is my comment
      val <- regexpr("\\s", name)
      if (val != -1) {
        name_comment <- substring(name, c(1, val + 1), c(val - 1, nchar(name)))
        return(name_comment[1])
      }
      return(name)
    },

    # create homopolymer of length size
    get_homop = function(size = 2) {
      as <- paste(rep("A", size), collapse = "")
      ts <- paste(rep("T", size), collapse = "")
      cs <- paste(rep("C", size), collapse = "")
      gs <- paste(rep("G", size), collapse = "")
      ns <- paste(rep("N", size), collapse = "")

      return(c(as, ts, cs, gs, ns))
    }
  )
)
