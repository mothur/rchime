#' @title sequence_dataset
#' @description 'sequence_dataset' is an R6 class that represents sequence FASTA
#'  and abundance data. It is a parent class for the 'sequence_data_table' and
#'   'SequenceDuckDB' classes whose data can be sourced from mothur formatted
#'    FASTA and count tables, XXXXX duckDB database XXXX, or created by
#'    mothur2's make_contigs function.
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
    #'   dataset <- sequence_data_table$new()
    #' @return A new `sequence_data_table` object which inherits
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
      num_non_ambig <- stri_count(seq, regex = "[*.*-ATGC]")
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
    #' Get location of first base in sequence
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_start = function(seq) {
      return(stringi::stri_locate_first_regex(seq, "[A-Z]")[1])
    },

    #' @description
    #' Get location of last base in sequence
    #' @param seq String, containing nucleotides
    #' @return An Integer
    calc_end = function(seq) {
      return(stringi::stri_locate_last_regex(seq, "[A-Z]")[1])
    },

    #' @description
    #' Read a mothur formatted count file
    #' @param names names of sequence in dataset (for matching sanity checks
    #'  between fasta and count files)
    #' @param filename Count file name (required)
    #' @param path path to count file (optional)
    #' @examples
    #'   dataset <- sequence_dataset$new()
    #'   fasta_data <- dataset$read_fasta_file(filename =
    #'     mothur2_example("test.fasta"))
    #'   fasta_names <- fasta_data$names
    #'   count_data <- dataset$read_count_file(names = fasta_names,
    #'    filename = mothur2_example("test.count_table"))
    #' @return A data.frame containing the sequence abundance data
    read_count_file = function(names, filename, path = NULL) {
      directory <- dirname(filename)
      filename <- basename(filename)

      # if path not set, make directory of filename
      if (is.null(path)) {
        path <- directory
      }

      count_file <- file.path(file.path(path), filename)
      file_conn <- file(count_file)
      file_data <- readLines(file_conn)
      num_lines <- length(file_data)
      close(file_conn)
      df <- data.frame()

      # check for group info in header
      if (num_lines > 1) {
        comment <- stri_locate_first_fixed(file_data[2], "#")

        # is this file compressed
        if (!is.na(comment[1])) {
          df <- private$fill_from_compressed_format(file_data, names)
        } else {
          df <- private$fill_from_uncompressed_format(file_data, names)
        }
      }
      return(df)
    },

    #' @description
    #' Read a FASTA formatted sequence file
    #' @param filename FASTA file name (required)
    #' @param path path to FASTA file (optional)
    #' @examples
    #'   dataset <- sequence_dataset$new()
    #'   fasta_data <- dataset$read_fasta_file(filename =
    #'   mothur2_example("test.fasta"))
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
          stri_count(x,
            regex = "[N]"
          )
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
          stri_locate_first_regex(
            x, "[A-Z]"
          )[1]
        })
      ))

      # calc ends
      ends <- unlist(lapply(
        df$Sequence,
        (function(x) {
          stri_locate_last_regex(
            x, "[A-Z]"
          )[1]
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
    #' Write count table data to mothur formatted count file
    #' @param df a data.table containing the count table data
    #' @param filename String, name of file you want to write to
    #' @param path path to count file (optional)
    #' @examples
    #'   dataset <- sequence_dataset$new()
    #'   fasta_data <- dataset$read_fasta_file(filename =
    #'   mothur2_example("test.fasta"))
    #'   fasta_names <- fasta_data$names
    #'   count_data <- dataset$read_count_file(names = fasta_names,
    #'    filename = mothur2_example("test.count_table"))
    #'   dataset$write_count_file(count_data,
    #'    filename = "output.count")
    #'   remove_file("output.count")
    write_count_file = function(df, filename, path = NULL) {
      if (file.exists(filename)) {
        file.remove(filename)
      }

      file_conn <- file(filename, open = "a")
      group_names <- names(df)[2:ncol(df)]
      has_group_data <- TRUE

      if (group_names[1] == "total") {
        group_names <- NULL
        has_group_data <- FALSE
      }
      writeLines(private$get_header(group_names), file_conn)

      ncols <- ncol(df)
      nrows <- nrow(df)
      if (has_group_data) {
        for (i in seq_len(nrows)) {
          # seq1.  10.  1,6 2,10 ....
          counts <- df[i, 2:ncols]
          counts <- c(sum(counts), counts)
          output <- paste(counts, collapse = "\t")
          output <- paste(df[i, 1], output, sep = "\t")
          writeLines(output, file_conn)
        }
      } else {
        for (i in seq_len(nrows)) {
          counts <- df[i, 2:ncols]
          total <- sum(counts)
          output <- paste(df[i, 1], total, sep = "\t")
          writeLines(output, file_conn)
        }
      }

      close(file_conn)
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
    #'   mothur2_example("test.fasta"))
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

    # Create header for mothur formatted count file
    get_header = function(groups = NULL) {
      if (!is.null(groups)) {
        output <- paste("Representative_Sequence     total   ")
        for (i in seq_along(groups)) {
          output <- paste(output, groups[i], sep = "\t")
        }
        return(output)
      } else {
        return("Representative_Sequence     total")
      }
    },

    # create homopolymer of length size
    get_homop = function(size = 2) {
      as <- paste(rep("A", size), collapse = "")
      ts <- paste(rep("T", size), collapse = "")
      cs <- paste(rep("C", size), collapse = "")
      gs <- paste(rep("G", size), collapse = "")
      ns <- paste(rep("N", size), collapse = "")

      return(c(as, ts, cs, gs, ns))
    },

    # Compressed Format: groupIndex,abundance. 1,6 means an ...
    # 	2,sample2	3,sample3	1,sample4
    # Representative_Sequence     total  sample2	sample3	sample4
    # returns data.frame contain count tables abundance data
    fill_from_compressed_format = function(file_data, fasta_names) {
      # extract group names
      # line 2 looks like: "#2,sample2	3,sample3	1,sample4"
      # remove first '#'
      comment <- stri_locate_first_fixed(file_data[2], "#")
      file_data[2] <- stri_sub(file_data[2], from = comment[1] + 1)
      words <- split_white_space(file_data[2])
      num_seqs <- length(fasta_names)

      has_group_data <- TRUE
      groups <- c()
      count_names <- rep("", num_seqs)
      df <- data.frame(names = count_names)

      for (i in seq_along(words)) {
        # parse group name
        file_index <- split_at_char(words[i], ",")

        # save group names
        groups <- c(groups, file_index[2])

        # prefill with 0
        samplei <- rep(0, num_seqs)
        df <- cbind(df, samplei)
      }
      names(df) <- c("names", groups)

      row <- 1
      # read compressed data lines
      for (i in 4:length(file_data)) {
        words <- split_white_space(file_data[i])

        name <- words[1]

        if (name %in% fasta_names) {
          # set name
          df[row, 1] <- name

          for (j in 3:length(words)) {
            data <- split_at_char(words[j], ",")
            # add all samples
            df[row, as.integer(data[1]) + 1] <- as.integer(data[2])
          }
          row <- row + 1
        } else {
          cli::cli_abort("[ERROR]: The sequence {.var {name}} is in
                    your count file and not in your dataset.")
        }
      }

      return(df)
    },

    # returns data.frame contain count tables abundance data
    fill_from_uncompressed_format = function(file_data, fasta_names) {
      # uncompressed format
      # Representative_Sequence  total  sample2	sample3	sample4
      words <- split_white_space(file_data[1])
      num_seqs <- length(fasta_names)
      has_group_data <- TRUE

      # no groups in file
      groups <- c()
      names <- rep("", num_seqs)
      df <- data.frame(names = names)

      if (length(words) == 2) {
        samplei <- rep(0, num_seqs)
        df <- cbind(df, samplei)
        groups <- c(groups, "total")
        has_group_data <- FALSE
      } else {
        for (i in 3:length(words)) {
          # save group names
          groups <- c(groups, words[i])

          # prefill with 0
          samplei <- rep(0, num_seqs)
          df <- cbind(df, samplei)
        }
      }
      names(df) <- c("names", groups)

      # read uncompressed data
      row <- 1
      for (i in 2:length(file_data)) {
        words <- split_white_space(file_data[i])

        if (length(words) >= 2) {
          name <- words[1]
          seq_total <- as.integer(words[2])

          if (name %in% fasta_names) {
            # set name
            df[row, 1] <- name

            if (has_group_data) {
              col <- 2
              for (j in 3:length(words)) {
                # add all samples
                df[row, col] <-
                  as.integer(words[j])
                col <- col + 1
              }
            } else {
              df[row, 2] <- as.integer(seq_total)
            }
            row <- row + 1
          } else {
            cli::cli_abort("[ERROR]: The sequence {.var {name}} is
                        in your count file and not in your dataset.")
          }
        }
      }
      return(df)
    }
  )
)
