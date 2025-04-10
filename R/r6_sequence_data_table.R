#' @title sequence_data_table
#' @description 'sequence_data_table' is an R6 class that represents sequence
#'  FASTA and abundance data. It inherits from 'sequence_dataset'.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom waldo compare
#' @import cli
#' @import r2r
#' @import data.table
#' @import stringi
#' @export
sequence_data_table <- R6Class("sequence_data_table",
  inherit = sequence_dataset,
  public = list(
    #' @description
    #' Create a new FASTA sequence dataset
    #' @param filename Name of file containing FASTA reads
    #' @param path Path to file containing FASTA reads
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'     rchime_example("test.fasta"))
    #'
    #' @return A new `sequence_data_table` object.
    initialize = function(filename = NULL, path = NULL) {
      super$initialize()
      self$clear()
      if (!is.null(filename)) {
        private$seq_data <- setDT(super$read_fasta_file(filename, path))

        num_table_rows <- 1
        all_names <- private$seq_data$names
        for (i in seq_along(all_names)) {
          private$seq2row[all_names[i]] <- num_table_rows
          num_table_rows <- num_table_rows + 1
        }

        # all seqs are assummed 'good'
        rows <- length(private$seq_data$names)
        private$table_seqs <- rep(TRUE, rows)

        abunds <- rep(1, rows)

        # add names with abundance 1
        private$add_seqs_abunds(names = private$seq_data$names, abunds = abunds)
        private$total <- nrow(private$seq_data)
      }
      invisible(self)
    },

    #' @description
    #' Get summary of sequences data - (summary_seqs output)
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$print()
    print = function() {
      self$print_bad_accnos()
      summary <- summary_seqs(dataset = self, silent = TRUE)
      print(summary)
      self$print_count_summary()
      cat(paste("\nNumber of unique seqs:", self$get_num_unique_seqs()), "\n")
      cat(paste("Total number of seqs:", self$get_num_seqs(), "\n"), "\n")
    },

    #' @description View summary of removed sequences
    #' 'print_bad_accnos()' outputs a summary of the reasons sequences were
    #'  removed from the dataset
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'     rchime_example("test.fasta"))
    #'
    #'   # flag 10 'bad' sequences
    #'   names <- dataset$get_names()[c(1:10)]
    #'   trash_codes <- rep("first_ten", 10)
    #'
    #'   dataset$remove_seqs(names, trash_codes)
    #'
    #'   dataset$print_bad_accnos()
    print_bad_accnos = function() {
      codes <- keys(private$accnos)
      cat("\nTrash_code   Unique_count    Total_count:\n")
      for (i in seq_along(codes)) {
        counts <- private$accnos[[codes[[i]]]]
        if (sum(counts) != 0) {
          cat(paste(codes[i], counts[1], counts[2], sep = "\t"), "\n")
        }
      }

      cat("\nNumber of unique removed seqs: ", private$unique_bad, sep = "")
      cat("\nTotal number of removed seqs: ", private$total_bad, "\n\n",
        sep = ""
      )
    },

    #' @description View summary of count table data
    #' 'print_count_summary()' outputs a summary of the sequence
    #'  abundance data
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'    rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   dataset$print_count_summary()
    print_count_summary = function() {
      if (private$has_group_data) {
        cat("Group   Total:\n")
        group_names <- self$get_groups()
        group_totals <- self$get_group_totals()
        for (i in seq_along(group_names)) {
          cat(paste(group_names[i], group_totals[i], sep = "\t"), "\n")
        }
      }
    },

    #' @description
    #' Add new sequence data
    #' @param names a vector of sequence names
    #' @param sequences a vector of sequence data
    #' @param comments a vector of sequence comments, (optional)
    #' @param olengths a vector containing sequences overlap lengths,
    #'  (optional)
    #' @param ostarts a vector containing sequences overlap starts, optional
    #' @param oends a vector containing sequences overlap ends, optional
    #' @param mismatches a vector containing the number of sequence
    #'  mismatches, (optional)
    #' @param ee a vector containing the sequences expected errors, optional
    #' @examples
    #'
    #'   dataset <- sequence_data_table$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   sequences <- c("ATGGGCT", "..TG--ACCGT..", "..GGuatgc..")
    #'   comments <- c("", "myComment", "")
    #'   dataset$add_seqs(names, sequences, comments)
    #'
    add_seqs = function(names, sequences, comments = c(),
                        olengths = c(), ostarts = c(),
                        oends = c(), mismatches = c(), ee = c()) {
      rows <- length(names)
      if (rows != 0) {
        # if no comments are created add them
        if (length(comments) == 0) {
          comments <- rep("", rows)
        }

        # innocent until proven guilty
        trash_codes <- rep("", rows)

        # calc numNs
        numns <- unlist(lapply(
          sequences,
          (function(x) stri_count(x, regex = "[N]"))
        ))
        # calc numBases
        lengths <- unlist(lapply(
          sequences,
          (function(x) super$calc_numbases(x))
        ))

        # calc starts
        starts <- unlist(lapply(sequences, (function(x) super$calc_start(x))))

        # calc ends
        ends <- unlist(lapply(sequences, (function(x) super$calc_end(x))))

        ambigs <- unlist(lapply(sequences, (function(x) super$calc_ambigs(x))))

        polymers <- unlist(lapply(
          sequences,
          (function(x) super$calc_longest_polymer(x))
        ))

        data <- data.table(
          names = names, sequences = sequences,
          comments = comments, trash_codes = trash_codes,
          starts = starts, ends = ends, lengths = lengths,
          ambigs = ambigs, polymers = polymers, numns = numns
        )

        # all rows and name, sequence and comments columns
        new_row <- data[c(1:rows), ]

        # add new rows to dataset and set new seqs status to "good"
        private$seq_data <- rbind(private$seq_data, new_row)
        private$table_seqs <- c(private$table_seqs, rep(TRUE, rows))

        num_table_rows <- 1
        all_names <- self$get_names()
        for (i in seq_along(all_names)) {
          private$seq2row[all_names[i]] <- num_table_rows
          num_table_rows <- num_table_rows + 1
        }

        # if no overlap_lengths are created add them
        if (length(olengths) != 0) {
          all_lengths <- c(
            length(olengths), length(ostarts), length(oends),
            length(mismatches), length(ee), length(names)
          )

          # if they are all the same length, then add
          if (length(unique(all_lengths)) == 1) {
            data <- data.table(
              olengths = olengths, ostarts = ostarts, oends = oends,
              mismatches = mismatches, ee = ee
            )

            # all rows and name, sequence and comments columns
            new_row <- data[c(1:rows), ]

            # add new rows to dataset
            private$contigs_data <- rbind(private$contigs_data, new_row)
          } else {
            cli::cli_alert("The contigs data length is mismatched, removing.")
          }
        }

        abunds <- rep(1, rows)

        # add names with abundance 1
        private$add_seqs_abunds(names = names, abunds = abunds)

        private$total <- sum(nrow(private$seq_data))
      }
      invisible(self)
    },

    #' @description
    #' Remove all sequences from dataset
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$clear()
    clear = function() {
      # fasta data
      private$seq_data <- data.table()
      # contigs data
      private$contigs_data <- data.table()
      # align data
      private$align_data <- data.table()
      # abundance data
      private$counts <- hashmap()
      # map to quick reference sequence data
      private$seq2row <- hashmap()
      # group names mapped to column in counts
      private$groups <- hashmap()
      # map of reasons for removal
      private$accnos <- hashmap()
      # map of pre-calculated group totals
      private$group_totals <- hashmap()

      # sequence counts
      private$total_bad <- 0
      private$unique_bad <- 0
      private$total <- 0

      # include groups and sequences
      private$table_groups <- c()
      private$table_seqs <- c()

      # is abundance data parsed by sample
      private$has_group_data <- FALSE

      invisible(self)
    },

    #' @description
    #' Get abundance for sequence in the dataset
    #' @param name String, Name of sequence
    #' @param group (optional) Name of group to get abundance for,
    #'   if not provided the total abundance for the sequence is returned
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   dataset$get_abund(name =
    #'    "M00967_43_000000000-A3JHG_1_1111_8697_7063")
    #' @return An integer
    get_abund = function(name, group = NULL) {
      seq_index <- private$seq2row[[name]]

      # if the name is NOT in the table
      if (is.null(seq_index)) {
        super$alert_missing_name(name, 0)
        return(0)
      }

      # is sequence 'active'
      if (private$table_seqs[seq_index]) {
        if (is.null(group)) {
          if (private$has_group_data) {
            return(sum(self$get_abunds(name)))
          } else {
            return(private$counts[[name]][[1]][1])
          }
        } else if (self$has_group(group)) {
          group_index <- private$groups[[group]]
          group_abund_data <- private$counts[[name]]
          seqs_inc_index <- which(group_abund_data[[1]] %in% group_index)
          if (length(seqs_inc_index) != 0) {
            return(group_abund_data[[2]][seqs_inc_index])
          }
        }
      } else {
        super$alert_eliminated_name(name, 0)
      }
      return(0)
    },

    #' @description
    #' Get abundances for sequence by group
    #' @param name String, Name of sequence
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   dataset$get_abunds(name =
    #'    "M00967_43_000000000-A3JHG_1_1111_8697_7063")
    #' @return Vector of abundances, in the same order as the groups
    #'             returned by dataset$get_groups().
    get_abunds = function(name) {
      seq_info <- self$in_table(name)
      if (seq_info$include) {
        if (private$has_group_data) {
          # c(3,2) c(100, 45)
          group_abund_data <- private$counts[[name]]
          abunds <- rep(0, length(private$groups))
          # group_abund_data[[1]] -> c(3,2)
          for (i in seq_along(group_abund_data[[1]])) {
            # group_abund_data[[2]] -> c(100,45)
            # i -> 1, abunds[3] <- 100
            # i -> 2, abunds[2] <- 45
            abunds[group_abund_data[[1]][i]] <- group_abund_data[[2]][i]
          }
          return(abunds[private$table_groups])
        } else {
          return(private$counts[[name]][[1]])
        }
      } else {
        super$alert_missing_name(name, c())
        return(c())
      }
    },

    #' @description
    #' Get data.table containing a report of sequences removed
    #' @examples
    #'  dataset <- sequence_data_table$new(filename =
    #'  rchime_example("test.fasta"))
    #'  some_seqs_to_remove <- c("M00967_43_000000000-A3JHG_1_2113_21675_17687",
    #'  "M00967_43_000000000-A3JHG_1_1114_17672_13068")
    #'  reasons <- c("example_removal", "example_removal")
    #'  dataset$remove_seqs(some_seqs_to_remove, reasons)
    #'  df <- dataset$get_accnos_report()
    get_accnos_report = function() {
      codes <- keys(private$accnos)
      unique_count <- c()
      total_count <- c()
      trash_codes <- c()
      for (i in seq_along(codes)) {
        counts <- private$accnos[[codes[[i]]]]

        if (sum(counts) != 0) {
          unique_count <- c(unique_count, counts[1])
          total_count <- c(total_count, counts[2])
          trash_codes <- c(trash_codes, codes[i])
        }
      }
      df <- data.table(
        trash_codes = trash_codes,
        unique_count = unique_count,
        total_count = total_count
      )

      return(list(
        summary = df, unique_bad = private$unique_bad,
        total_bad = private$total_bad
      ))
    },

    #' @description
    #' Get data.table containing align report data used for screening
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'     rchime_example("test.fasta"))
    #'   reference <- sequence_data_table$new(filename =
    #'     rchime_example("reference.fasta"))
    #'
    #'   # results$align_report contains complete align report
    #'   # results$accnos will contain the names of any sequences that were
    #'   #   flipped
    #'
    #'   # results <- align_seqs(dataset, reference)
    #'
    #'   # report contains partial align report used for screening:
    #'   # search_score, sim_score and longest_insert
    #'   report <- dataset$get_align_report()
    #'
    get_align_report = function() {
      if (self$has_align_report()) {
        searches <- private$align_data$search_scores[private$table_seqs]
        sims <- private$align_data$sim_scores[private$table_seqs]
        inserts <- private$align_data$longest_inserts[private$table_seqs]

        align_report <- data.table(
          search_scores = searches,
          sim_scores = sims,
          longest_inserts = inserts
        )

        return(align_report)
      }
      return(data.table())
    },

    #' @description
    #' Get data.table containing contigs report data
    #' @examples
    #'   # You can use the mothur2 package to assemble your
    #'   # paired reads and create a dataset with contigs data
    #'
    #'   # files <- make_file(input_dir = rchime_example(), type = "gz")
    #'   # contigs <- make_contigs(files)
    #'
    #'   # df <- contigs$dataset$get_contigs_report()
    #'
    get_contigs_report = function() {
      if (self$has_contigs_report()) {
        lengths <- private$seq_data$lengths[private$table_seqs]
        olengths <- private$contigs_data$olengths[private$table_seqs]
        ostarts <- private$contigs_data$ostarts[private$table_seqs]
        oends <- private$contigs_data$oends[private$table_seqs]
        mismatches <- private$contigs_data$mismatches[private$table_seqs]
        numns <- private$seq_data$numns[private$table_seqs]
        ee <- private$contigs_data$ee[private$table_seqs]

        contigs_report <- data.table(
          lengths = lengths,
          olengths = olengths,
          ostarts = ostarts,
          oends = oends,
          mismatches = mismatches,
          numns = numns,
          ee = ee
        )

        return(contigs_report)
      }
      return(data.table())
    },

    #' @description
    #' Get count table data
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   count_table <- dataset$get_count_table()
    #' @return a data.table containing the count table data
    get_count_table = function() {
      df <- data.table(names = self$get_names())
      groups <- c()
      if (private$has_group_data) {
        # add group columns
        groups <- self$get_groups()
        for (i in seq_along(groups)) {
          # fill with 0
          samplei <- rep(0, length(names))
          df <- cbind(df, samplei)
        }

        names <- self$get_names()
        for (i in seq_along(names)) {
          counts <- self$get_abunds(names[i])
          for (j in seq_along(counts)) {
            set(df, i, as.integer(j + 1), counts[j])
          }
        }
        names(df) <- c("names", groups)
      } else {
        df <- cbind(df, self$get_seqs_abunds())
        names(df) <- c("names", "total")
      }
      return(df)
    },

    #' @description
    #' Get names of groups in the dataset
    #' @param name The name of the sequence you want groups for, optional
    #' @return A character vector
    get_groups = function(name = NULL) {
      if (is.null(name)) {
        return(private$get_included_groups())
      } else {
        seq_info <- self$in_table(name)
        if (seq_info$include && private$has_group_data) {
          inc_groups <- private$get_included_groups()
          group_abund_data <- private$counts[[name]]
          included_indexes <- private$get_included_groups_indexes()
          seqs_group_inds <- which(group_abund_data[[1]] %in% included_indexes)
          # returns groups that are included in dataset and present in seq
          return(inc_groups[group_abund_data[[1]][seqs_group_inds]])
        }
        return(c())
      }
    },

    #' @description
    #' Get the number of sequences represented in samples
    #' @param group String, name of sample
    #' @return A vector of integers
    get_group_totals = function(group = NULL) {
      if (!is.null(group)) {
        group_total <- private$group_totals[[group]]
        if (!is.null(group_total)) {
          return(group_total)
        }
      } else if (private$has_group_data) {
        groups <- self$get_groups()
        totals <- rep(0, length(groups))

        for (i in seq_along(groups)) {
          totals[i] <- private$group_totals[[groups[i]]]
        }

        return(totals)
      }
      return(0)
    },

    #' @description
    #' Get names of sequences in the dataset
    #' @param group String, name of sample
    #' @examples
    #'   dataset <- sequence_data_table$new(
    #'   filename = rchime_example('test.fasta'))
    #'
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'
    #'   # get names of sequences in your dataset
    #'   names <- dataset$get_names()
    #'
    #'   # get names of sequences in sample 'F3D0'
    #'   names <- dataset$get_names(group = 'F3D0')
    #'
    get_names = function(group = NULL) {
      # invalid groups will produce null index
      group_index <- private$groups[[group]]

      if (!is.null(group_index)) {
        return(private$select_group(
          private$seq_data$names,
          group_index
        ))
      } else if (!is.null(group)) {
        cli::cli_alert("{.var {group}} is invalid, ignoring.")
      }
      return(private$seq_data$names[private$table_seqs])
    },

    #' @description
    #' Get number of groups in the dataset
    #' @return An integer
    get_num_groups = function() {
      if (private$has_group_data) {
        return(sum(private$table_groups, na.rm = TRUE))
      }
      return(0)
    },

    #' @description
    #' Get the number of sequences in the dataset
    #' @param group The name of the group you want number of sequences for,
    #'  optional
    #' @return An integer
    get_num_seqs = function(group = NULL) {
      if (is.null(group)) {
        return(private$total)
      } else {
        return(self$get_group_totals(group))
      }
      return(0)
    },

    #' @description
    #' Get number of unique sequences in dataset
    #' @return An integer
    get_num_unique_seqs = function() {
      return(sum(private$table_seqs, na.rm = TRUE))
    },

    #' @description
    #' Get data.table containing FASTA data
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   seqs <- dataset$get_seqs_table()
    #' @return data.table
    get_seqs_table = function() {
      return(private$get_good())
    },

    #' @description
    #' Get vector containing FASTA nucleotide strings
    #' @param group String, name of sample
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'
    #'   # to get all sequences in dataset
    #'   seqs <- dataset$get_seqs()
    #'
    #'   # to get sequences from sample 'F3D0'
    #'   seqs <- dataset$get_seqs(group = 'F3D0')
    #'
    #' @return data.table
    get_seqs = function(group = NULL) {
      # invalid groups will produce null index
      group_index <- private$groups[[group]]

      if (!is.null(group_index)) {
        return(private$select_group(
          private$seq_data$sequences,
          group_index
        ))
      } else if (!is.null(group)) {
        cli::cli_alert("{.var {group}} is invalid, ignoring.")
      }
      return(private$seq_data$sequences[private$table_seqs])
    },

    #' @description
    #' Get abundances for sequences in dataset
    #' @param group String, name of sample
    #' @examples
    #'   dataset <- sequence_data_table$new(
    #'   filename = rchime_example('test.fasta'))
    #'
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'
    #'   # get abundances of all sequences in your dataset
    #'   abunds <- dataset$get_seqs_abunds()
    #'
    #'   # get abundances of sequences in sample 'F3D0'
    #'   abunds <- dataset$get_seqs_abunds(group = 'F3D0')

    #' @return Vector of abundances, in the same order as the sequences
    #'             returned by dataset$get_names().
    get_seqs_abunds = function(group = NULL) {
      abunds <- c()
      if (private$has_count()) {
        # gave us an invalid group
        if (!is.null(group)) {
          if (is.null(private$groups[[group]])) {
            cli::cli_alert("{.var {group}} is invalid, ignoring.")
            group <- NULL
          }
        }
        names <- self$get_names(group)
        for (i in seq_along(names)) {
          abunds[i] <- self$get_abund(names[i], group)
        }
      } else {
        abunds <- rep(1, self$get_num_unique_seqs())
      }

      return(abunds)
    },

    #' @description
    #' Get data.table containing summary report data
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   df <- dataset$get_summary_report()
    get_summary_report = function() {
      lengths <- private$seq_data$lengths[private$table_seqs]
      ambigs <- private$seq_data$ambigs[private$table_seqs]
      polymers <- private$seq_data$polymers[private$table_seqs]
      numns <- private$seq_data$numns[private$table_seqs]
      starts <- private$seq_data$starts[private$table_seqs]
      ends <- private$seq_data$ends[private$table_seqs]

      summary_report <- data.table(
        lengths = lengths,
        starts = starts,
        ends = ends,
        ambigs = ambigs,
        polymers = polymers,
        numns = numns
      )

      return(summary_report)
    },

    #' @description
    #' Determine if the dataset contains align report data
    #' @return Boolean
    has_align_report = function() {
      if (nrow(private$align_data) == 0) {
        return(FALSE)
      }
      return(TRUE)
    },

    #' @description
    #' Determine if the dataset contains contigs report data
    #' @return Boolean
    has_contigs_report = function() {
      if (nrow(private$contigs_data) == 0) {
        return(FALSE)
      }
      return(TRUE)
    },

    #' @description
    #' Determine if a group is present in the dataset
    #' @param group String, Name of sample
    #' @return Boolean
    has_group = function(group) {
      if (private$has_group_data) {
        group_index <- private$groups[[group]]

        # in case the group is not in table
        if (!is.null(group_index)) {
          return(private$table_groups[group_index])
        }
      }
      return(FALSE)
    },

    #' @description
    #' Determine if a sequence is present in the dataset and its index
    #' @param name String, Name of sequence
    #' @return list, list(include = boolean, index = row in table)
    in_table = function(name) {
      seq_index <- private$seq2row[[name]]
      # if the name is in the table
      if (!is.null(seq_index)) {
        return(list(include = TRUE, index = seq_index))
      }
      return(list(include = FALSE, index = seq_index))
    },

    #' @description
    #' Determine if the dataset is aligned and its alignment length
    #' @return list(aligned = (TRUE OR FALSE),
    #'              length = (alignment length or NA))
    is_aligned = function() {
      results <- list(aligned = FALSE, length = NA)
      seq_lengths <- private$get_seqs_lengths()

      if (length(unique(seq_lengths)) == 1) {
        results$aligned <- TRUE
        results$length <- seq_lengths[1]
      } else {
        results$length <- unique(seq_lengths)
      }

      return(results)
    },

    #' @description merges identical sequences in the dataset
    #' @param seqs_to_merge List of vectors, each vector will be merged
    #' @param group String, when a group is provided only that groups counts
    #'                       will be merged (optional)
    #' @param reason String, reason for merging sequence
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'
    #'   # Number of Unique seqs: 1000
    #'   # Total number of seqs: 1000
    #'
    #'   seqs_to_merge <- list(c("M00967_43_000000000-A3JHG_1_1112_15471_15261",
    #'   "M00967_43_000000000-A3JHG_1_2114_21028_11825",
    #'   "M00967_43_000000000-A3JHG_1_2106_5820_9300"),
    #'   c("M00967_43_000000000-A3JHG_1_1111_26505_10426",
    #'   "M00967_43_000000000-A3JHG_1_1111_23613_8947"))
    #'
    #'   dataset$merge_seqs(seqs_to_merge)
    #'
    #'   # Number of Unique seqs: 997
    #'   # Total number of seqs: 1000
    merge_seqs = function(seqs_to_merge, reason = "merged", group = NULL) {
      lapply(seqs_to_merge, (function(x) {
        private$merge_seq(x,
          reason = reason,
          group = group
        )
      }))
      invisible(self)
    },

    #' @description
    #' Reinstates sequences removed for the reasons in trash_codes
    #' @param trash_codes vector containing reasons for sequences removal
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'
    #'   # flag 10 sequences for removal
    #'   names <- dataset$get_names()[c(1:10)]
    #'   trash_codes <- rep("first_ten", 10)
    #'
    #'   dataset$remove_seqs(names, trash_codes)
    #'
    #'   # reinclude sequences you removed
    #'   dataset$reinstate_seqs(trash_codes = c("first_ten"))
    reinstate_seqs = function(trash_codes) {
      if (length(trash_codes) != 0) {
        # for each trash code
        for (i in seq_along(trash_codes)) {
          trash_code <- trash_codes[i]
          # find sequences with that trash code
          seqs_with_trashcode <- stri_detect_fixed(
            private$seq_data$trash_codes,
            trash_code
          )

          # not all seqs with this trash code will be reinstated. The
          # sequence may have more than one trashcode
          names_2add <- private$get_seqs_for_reinstate(trash_code)

          accnos_names <- private$seq_data$names[seqs_with_trashcode]

          groups <- self$get_groups()

          groups_to_re_add <- c()

          for (j in seq_along(names_2add)) {
            seq_row <- private$seq2row[[names_2add[j]]]

            # set include seq flags
            private$table_seqs[seq_row] <- TRUE
            private$seq_data$trash_codes[seq_row] <- ""

            # check if this seq re-adds a group
            this_seqs_groups <- private$get_groups_for_reinstate(names_2add[j])

            groups_to_re_add <- setdiff(this_seqs_groups, groups)

            if (length(groups_to_re_add) != 0) {
              for (k in seq_along(groups_to_re_add)) {
                group_index <- private$groups[[groups_to_re_add[k]]]
                private$table_groups[group_index] <- TRUE
              }
              groups <- self$get_groups()
            }

            abunds <- self$get_abunds(names_2add[j])
            if (private$has_group_data) {
              diff_abunds <- abunds * -1L
              private$update_group_totals(diff_abunds)
            }
            abund <- sum(abunds)

            # update totals
            private$total <- private$total + abund
            private$unique_bad <- private$unique_bad - 1
            private$total_bad <- private$total_bad - abund
          }

          # update accnos counts
          for (j in seq_along(accnos_names)) {
            accnos_counts <- c(1, self$get_abund(accnos_names[j]))

            new_total <- private$accnos[[trash_code]] - accnos_counts
            if (sum(new_total) != 0) {
              private$accnos[[trash_code]] <- new_total
            } else {
              delete(private$accnos, trash_code)
            }
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Remove sequences from dataset for cause
    #' @param names vector containing names of sequences to remove
    #' @param trash_codes vector containing reasons for sequences removal
    #' @examples
    #'   dataset <- sequence_data_table$new(filename =
    #'                     rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'                     rchime_example("test.count_table"))
    #'   dataset
    #'
    #' # remove first 10 sequences from the dataset
    #'
    #'   names <- dataset$get_names()[c(1:10)]
    #'   trash_codes <- rep("first_ten", 10)
    #'
    #'   dataset$remove_seqs(names, trash_codes)
    #'   dataset
    #'
    remove_seqs = function(names, trash_codes) {
      if (length(names) != length(trash_codes)) {
        super$abort_length_mismatch(
          "names", "trash_codes", length(names),
          length(trash_codes)
        )
      } else {
        if (length(names) != 0) {
          groups <- self$get_groups()

          for (i in seq_along(names)) {
            row_to_remove <- private$seq2row[[names[i]]]
            if (!is.null(row_to_remove)) {
              private$remove_seq(names[i], trash_codes[i], TRUE)
            }
          }
        }

        # are their samples to remove now? This can happen if the seqs
        # removed eliminate a group
        new_totals <- self$get_group_totals()

        # are there groups in the dataset with no sequences
        if (0 %in% new_totals) {
          groups <- self$get_groups()
          for (i in seq_along(groups)) {
            if (new_totals[i] == 0) {
              group_index <- private$groups[[groups[i]]]
              private$table_groups[group_index] <- FALSE
            }
          }
        }
      }
      invisible(self)
    },
    #' @description
    #' Add group assignment data
    #' @param names a vector of sequence names
    #' @param groups a vector of group assignments
    #' @param filename String, Name of mothur formatted count file
    #' @param path String, Path to mothur formatted count file
    #' @examples
    #'
    #'   dataset <- sequence_data_table$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   sequences <- c("ATGGGCT", "..TG--ACCGT..", "..GGuatgc..")
    #'   groups <- c("sample1", "sample2", "sample3")
    #'   dataset$add_seqs(names, sequences)
    #'   dataset$set_group_assignments(names, groups)
    #'
    set_group_assignments = function(names = NULL, groups = NULL,
                                     filename = NULL, path = NULL) {
      # clear old assignments, if any
      private$counts <- hashmap()
      private$table_groups <- c()
      private$has_group_data <- FALSE
      private$group_totals <- hashmap()
      private$groups <- hashmap()

      if (!is.null(names) && !is.null(groups)) {
        if (length(names) != length(groups)) {
          super$abort_length_mismatch(
            "names", "groups", length(names),
            length(groups)
          )
        }
        if (length(names) != length(private$seq_data$names)) {
          super$abort_length_mismatch(
            "names", "groups", length(names),
            length(private$seq_data$names)
          )
        }

        for (i in seq_along(names)) {
          seq_info <- self$in_table(names[i])
          if (!seq_info$include) {
            super$abort_name_missing(names[i])
          }
        }

        # add groups
        unique_groups <- unique(groups)

        # are all the sequence assigned to "NA"
        if (length(unique_groups) > 1) {
          private$add_groups(unique_groups)
        } else if ((length(unique_groups) == 1) && unique_groups[1] != "NA") {
          private$add_groups(unique_groups)
        } else {
          private$has_group_data <- FALSE
        }

        num_groups <- length(unique_groups)

        if (num_groups != 0) {
          for (i in seq_along(groups)) {
            group_index <- private$groups[[groups[i]]]
            private$counts[[names[i]]] <- list(c(group_index), c(1))
          }
        }
      } else if (!is.null(filename)) {
        df <- super$read_count_file(self$get_names(), filename, path)
        group_names <- names(df)[2:ncol(df)]
        private$has_group_data <- TRUE

        if (nrow(df) != self$get_num_unique_seqs()) {
          super$abort_length_mismatch(
            "names", "count table rows", length(self$get_num_unique_seqs()),
            nrow(df)
          )
        }

        if (group_names[1] == "total") {
          group_names <- NULL
          private$has_group_data <- FALSE
          # add flag for "total"
          private$table_groups <- c(TRUE)
        } else {
          private$add_groups(group_names)
        }

        ncols <- ncol(df)
        for (i in seq_len(nrow(df))) {
          counts <- as.numeric(df[i, 2:ncols])
          group_indexes <- which(counts != 0)
          name <- as.character(df[i, 1])
          if (private$has_group_data) {
            private$counts[[name]] <-
              list(group_indexes, counts[group_indexes])
          } else {
            private$counts[[name]] <- list(counts[group_indexes])
          }
        }
      }

      # this count file has groups
      if (private$has_group_data) {
        sums <- rep(0, self$get_num_groups())
        names <- self$get_names()
        for (i in seq_along(names)) {
          counts <- self$get_abunds(names[i])
          for (i in seq_along(counts)) {
            sums[i] <- sums[i] + counts[i]
          }
        }

        # include all groups
        groups <- self$get_groups()
        private$table_groups <- c(rep(TRUE, length(groups)))

        for (i in seq_along(groups)) {
          private$group_totals[groups[i]] <- sums[i]
        }
      }

      private$calc_total()

      invisible(self)
    },

    #' @description
    #' Set abundance data for sequences
    #' @param names a vector of sequence names
    #' @param abunds a list of abundance assignments
    #' @param reason, String (optional) reason for setting sequences to
    #'  zero abundance
    #' @examples
    #'   dataset <- sequence_data_table$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   sequences <- c("ATGGGCT", "..TG--ACCGT..", "..GGuatgc..")
    #'
    #'   dataset$add_seqs(names, sequences)
    #'
    #'   # set total abundance. seq1 -> 100, seq2 -> 50, seq3 -> 5
    #'   abunds <- list(100, 50, 5)
    #'   dataset$set_abundances(names, abunds)
    #'
    #'   # add group data
    #'   groups <- c("sample1", "sample1", "sample2")
    #'   dataset$set_group_assignments(names, groups)
    #'
    #'   # set the abundance of 'seq1' in 'sample1' to 150
    #'   dataset$set_abundances(c("seq1"), list(c(150, 0)))
    #'
    #'   # set the abundance of 'seq1' in 'sample1' to 150,
    #'   #                      'seq1' in 'sample2' to 0 and
    #'   # set the abundance of 'seq3' in 'sample1' to 75,
    #'   #                      'seq3' in 'sample2' to 5
    #'   dataset$set_abundances(c("seq1", "seq3"), list(c(150, 0), c(75, 5)))
    #'
    set_abundances = function(names = NULL, abunds = NULL, reason = "merged") {
      if (!is.null(names) && !is.null(abunds)) {
        if (length(names) != length(abunds)) {
          super$abort_length_mismatch(
            "names", "abunds", length(names),
            length(abunds)
          )
        }

        reason <- paste(reason, ",", sep = "")

        for (i in seq_along(names)) {
          if (!self$in_table(names[i])$include) {
            super$abort_name_missing(names[i])
          } else {
            # self$get_abunds(names[i])
            if (sum(abunds[[i]]) == 0) {
              private$remove_seq(names[i], reason, TRUE)
            } else {
              # save old abunds if you have groups so we can update group_totals
              orig_abunds <- c()
              if (private$has_group_data) {
                orig_abunds <- self$get_abunds(names[i])
              }

              # update counts
              private$counts[[names[i]]] <- private$convert_to_sparse(abunds[i])

              if (private$has_group_data) {
                diff_abunds <- orig_abunds - abunds[[i]]

                # update group_totals
                private$update_group_totals(diff_abunds)
              }
            }
          }
        }

        if (private$has_group_data) {
          # remove groups that are now zeroed out
          new_totals <- self$get_group_totals()

          # are there groups in the dataset with no sequences
          if (0 %in% new_totals) {
            groups <- self$get_groups()

            for (i in seq_along(groups)) {
              if (new_totals[i] == 0) {
                group_index <- private$groups[[groups[i]]]
                private$table_groups[group_index] <- FALSE
              }
            }
          }
        }

        private$calc_total()
      } else {
        cli::cli_abort("You must provide both names and abunds.")
      }

      invisible(self)
    },

    #' @description
    #' Set existing sequence data
    #' @param names a vector of sequence names
    #' @param sequences a vector of sequence data
    #' @param comments a vector of sequence comments, (optional)
    #' @param search_scores a vector of alignment search scores,
    #'  (optional, used by align_seqs)
    #' @param sim_scores a vector of alignment similarity to template scores,
    #'  (optional, used by align_seqs)
    #' @param longest_inserts a vector of longest inserts in the alignment,
    #'  (optional, used by align_seqs)
    #'
    #' @examples
    #'   dataset <- sequence_data_table$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   sequences <- c("ATGGGCT", "..TG--ACCGT..", "..GGuatgc..")
    #'   dataset$add_seqs(names, sequences)
    #'
    #'   # change sequence string for seq2 from ..TG--ACCGT.. to ...ATG-G-GCT..
    #'   dataset$set_seqs(c("seq2"), c("...ATG-G-GCT.."))
    #'
    set_seqs = function(names, sequences, comments = c(),
                        search_scores = c(), sim_scores = c(),
                        longest_inserts = c()) {
      rows <- length(sequences)

      if (rows != 0) {
        # calc numNs
        numns <- unlist(lapply(
          sequences,
          (function(x) stri_count(x, regex = "[N]"))
        ))

        # calc numBases
        lengths <- unlist(lapply(
          sequences,
          (function(x) super$calc_numbases(x))
        ))

        # calc starts
        starts <- unlist(lapply(
          sequences,
          (function(x) super$calc_start(x))
        ))

        # calc ends
        ends <- unlist(lapply(sequences, (function(x) super$calc_end(x))))

        ambigs <- unlist(lapply(
          sequences,
          (function(x) super$calc_ambigs(x))
        ))

        polymers <- unlist(lapply(
          sequences,
          (function(x) super$calc_longest_polymer(x))
        ))

        cli_progress_bar("Updating sequences",
          total = length(names),
          clear = TRUE
        )

        indexes <- rep(-1, length(names))
        for (i in seq_along(names)) {
          indexes[i] <- private$seq2row[[names[i]]]

          cli_progress_update()

          private$seq_data[indexes[i], 2] <- sequences[i]
          private$seq_data[indexes[i], 5] <- starts[i]
          private$seq_data[indexes[i], 6] <- ends[i]
          private$seq_data[indexes[i], 7] <- lengths[i]
          private$seq_data[indexes[i], 8] <- ambigs[i]
          private$seq_data[indexes[i], 9] <- polymers[i]
          private$seq_data[indexes[i], 10] <- numns[i]

          if (length(comments) != 0) {
            private$seq_data[indexes[i], 3] <- comments[i]
          }
        }
        cli_progress_done()


        # if no overlap_lengths are created add them
        if (length(search_scores) != 0) {
          all_lengths <- c(
            length(search_scores), length(sim_scores),
            length(longest_inserts), length(names)
          )

          # if they are all the same length, then create align_data
          if (length(unique(all_lengths)) == 1) {
            cli_progress_bar("Adding align report data ",
              total = length(names),
              clear = TRUE
            )
            # create blank table, if needed
            if (!self$has_align_report()) {
              num_all_seqs <- length(private$table_seqs)

              data <- data.table(
                search_scores = rep(0, num_all_seqs),
                sim_scores = rep(0, num_all_seqs),
                longest_inserts = rep(0, num_all_seqs)
              )

              # all rows and name, sequence and comments columns
              new_row <- data[c(1:num_all_seqs), ]

              # add new rows to dataset
              private$align_data <- rbind(private$align_data, new_row)
            }

            # fill in "good" seqs align data
            for (i in seq_along(indexes)) {
              cli_progress_update()

              private$align_data[indexes[i], 1] <- search_scores[i]
              private$align_data[indexes[i], 2] <- sim_scores[i]
              private$align_data[indexes[i], 3] <- longest_inserts[i]
            }
            cli_progress_done()
          } else {
            cli::cli_alert("The align data length is mismatched, removing.")
          }
        }
      }
      invisible(self)
    }
  ),
  private = list(

    # FASTA data
    # names = c(), sequences = c(), comments = c(),
    # trash_codes = c(), starts = c(), ends = c(),
    # lengths = c(), ambigs = c(), polymers = c(), numns = c(),
    seq_data = data.table(),

    # contigs report data
    # olengths = c(), ostarts = c(), oends = c(),
    # mismatches = c(), ee = c() same order as seq_data
    contigs_data = data.table(),

    # align report data
    # search_score, sim_score and longest_insert
    align_data = data.table(),

    # count table data - sparse
    # sequence name -> list of group indexes and abundances.
    # data.table -> groups is a vector of integers representing the
    # .                            groups the sequence is included in
    #               abunds is a vector of integers representing the abundances
    # .                        of the groups the sequence is included in
    # groups['F3D2'] -> 1, groups['F3D78'] -> 2, groups['F3DN'] -> 3
    # groups['F3D1'] -> 4,  groups['F3D6'] -> 5,  groups['F3D12'] -> 6
    # seq1 -> list(c(4, 2, 6), c(2, 5, 3)) seq1 has a total abundance of 10,
    # 2 from 'F3D1', 5 from 'F3D78', 3 from 'F3D1'
    counts = hashmap(),

    # groupName -> index in private$table_groups
    groups = hashmap(),

    # sequence name -> row in seq_data
    # seq2row["seq1"] will return row in seq_data where seq1's info is
    seq2row = hashmap(),

    # vector of booleans indicating whether sample is present in 'current'
    # table.
    table_groups = c(),

    # vector of booleans indicating whether sequence is present in 'current'
    # table.
    table_seqs = c(),

    # trash_code -> count of seqs
    accnos = hashmap(),
    total_bad = 0,
    unique_bad = 0,
    total = 0,

    # is abundance data parsed by sample
    has_group_data = FALSE,

    # groupName -> total abundance for group
    # sample4 -> 250 means sample4 contains 250 sequences
    group_totals = hashmap(),
    add_seqs_abunds = function(names, abunds) {
      if (private$has_group_data) {
        cli::cli_abort("[ERROR]: The 'sequence_data_table' contains group
                      information, you must provide by sample data.")
      } else {
        for (i in seq_along(names)) {
          private$counts[[names[i]]] <- list(abunds[i])
        }
        private$table_groups <- c(TRUE)
      }
      invisible(self)
    },
    add_groups = function(groups) {
      groups <- sort(groups)
      private$table_groups <- c(rep(TRUE, length(groups)))

      for (i in seq_along(groups)) {
        private$groups[groups[i]] <- i
      }
      private$has_group_data <- TRUE

      invisible(self)
    },
    calc_total = function() {
      # dataset with abundance data parsed by sample
      if (private$has_group_data) {
        private$total <- sum(self$get_group_totals())
      } else if (private$has_count()) {
        # dataset with abundance data, no samples
        private$total <- sum(self$get_seqs_abunds())
      } else {
        # dataset without abundance data
        private$total <- self$get_num_unique_seqs()
      }

      invisible(self)
    },

    # abunds a vector of abundances in full format, convert to sparse format
    # seq1 -> c(0, 5, 0, 2, 0, 3) becomes
    # seq1 -> list(c(4, 2, 6), c(2, 5, 3))
    # 2 from 'sample4', 5 from 'sample2', 3 from 'sample6'
    convert_to_sparse = function(abunds) {
      if (private$has_group_data) {
        if (length(abunds[[1]]) != self$get_num_groups()) {
          cli::cli_abort("When setting abunds on a dataset with groups you
                             must provide an abundance for each sample.")
        }
        indexes <- which(abunds[[1]] != 0)
        groups <- self$get_groups()[indexes]
        group_indexes <- c()
        for (i in seq_along(groups)) {
          group_indexes <- c(group_indexes, private$groups[[groups[i]]])
        }

        return(list(group_indexes, abunds[[1]][indexes]))
      }
      return(abunds)
    },

    # Clear sequences from dataset
    finalize = function() {
      self$clear()
    },

    # Get names of groups in the dataset
    get_included_groups = function() {
      return(sort(unlist(keys(private$groups))[private$table_groups]))
    },
    get_included_groups_indexes = function() {
      return(which(private$table_groups))
    },
    get_good = function() {
      return(private$seq_data[private$seq_data$trash_codes == "", ])
    },
    get_groups_for_reinstate = function(name) {
      all_groups <- sort(unlist(keys(private$groups)))
      indexes <- private$counts[[name]][[1]]
      return(all_groups[indexes])
    },
    get_seqs_lengths = function() {
      return(unlist(lapply(
        private$seq_data$sequences,
        (function(x) stri_length(x))
      )))
    },

    # returns vector of names to add back into active state
    get_seqs_for_reinstate = function(trash_code) {
      # remove trash code from that sequences list of trash codes
      private$seq_data$trash_codes <-
        stri_replace_all_fixed(
          private$seq_data$trash_codes,
          trash_code, ""
        )

      # remove double comma
      private$seq_data$trash_codes <-
        stri_replace_all_fixed(private$seq_data$trash_codes, ",,", ",")

      # if it's now good, add the seq back in
      seqs_to_re_add <- (private$seq_data$trash_codes == ",")

      # get the names
      reinstated_names <- private$seq_data$names[seqs_to_re_add]
    },

    # Determine if the dataset contains abundance data
    # return Boolean
    has_count = function() {
      if (length(private$counts) == 0) {
        return(FALSE)
      }
      return(TRUE)
    },
    merge_seq = function(names, reason, group = NULL) {
      # make sure seqs are in the dataset
      for (j in seq_along(names)) {
        seq_index <- private$seq2row[[names[j]]]
        if (is.null(seq_index)) {
          super$abort_name_missing(names[j])
        }
      }

      group_index <- private$groups[[group]]
      reason <- paste(reason, ",", sep = "")

      if (length(names) != 1) {
        keeper_name <- names[1]
        names <- names[-c(1)]
        k_abunds <- self$get_abunds(keeper_name)

        for (i in seq_along(names)) {
          dup_name <- names[i]
          dup_abunds <- self$get_abunds(dup_name)

          if (is.null(group_index)) {
            # merge all counts
            k_abunds <- k_abunds + dup_abunds
            private$remove_seq(dup_name, reason, FALSE)
          } else {
            # merge counts for this group
            k_abunds[group_index] <- k_abunds[group_index] +
              dup_abunds[group_index]
            dup_abunds[group_index] <- 0
            if (sum(dup_abunds) == 0) {
              private$remove_seq(dup_name, "duplicate,", FALSE)
            } else {
              d <- list(dup_abunds)
              private$counts[[dup_name]] <- private$convert_to_sparse(d)
            }
          }
        }

        # add dups groups and abunds to keeper
        k_abunds <- list(k_abunds)
        private$counts[[keeper_name]] <- private$convert_to_sparse(k_abunds)
      }

      invisible(self)
    },
    remove_seq = function(name, trash_code, update_totals = TRUE) {
      # when merging sequences we don't want to update group totals
      abunds <- self$get_abunds(name)
      abund <- sum(abunds)

      # we are merging seqs, so don't update group totals, or numSeqs
      if (update_totals) {
        # remove seq from group totals
        if (private$has_group_data) {
          private$update_group_totals(abunds)
        }
        private$total <- private$total - abund
      }

      # add reason for removal
      these_reasons <- split_at_char(trash_code, ",")

      for (k in seq_along(these_reasons)) {
        if (has_key(private$accnos, these_reasons[k])) {
          this_seqs_counts <- c(1, abund)

          private$accnos[[these_reasons[k]]] <-
            private$accnos[[these_reasons[k]]] + this_seqs_counts
        } else {
          private$accnos[[these_reasons[k]]] <- c(1, abund)
        }
      }

      # add trash_code reasons to accnos totals
      private$unique_bad <- private$unique_bad + 1
      private$total_bad <- private$total_bad + abund

      # trash_code column = 4
      row_to_remove <- private$seq2row[[name]]
      private$seq_data[row_to_remove, 4] <- trash_code
      private$table_seqs[row_to_remove] <- FALSE

      invisible(self)
    },

    # passed seq names or sequences and returns items in table and in group
    select_group = function(table_column, group_index) {
      column_entries_in_group <- c()

      if (private$table_groups[group_index]) {
        all_entries <- table_column[private$table_seqs]
        names <- private$seq_data$names[private$table_seqs]

        for (i in seq_along(names)) {
          group_abunds <- private$counts[[names[i]]]
          if (group_index %in% group_abunds[[1]]) {
            column_entries_in_group <- c(
              column_entries_in_group,
              all_entries[i]
            )
          }
        }
      }
      return(column_entries_in_group)
    },

    # diff_abunds -> difference between old abunds and new abunds
    # seq_old_abunds c(20, 10, 40)
    # seq_new_abunds c(0, 20, 35)
    # diff_abunds.   c(20, -10, 5)
    update_group_totals = function(diff_abunds) {
      if (sum(diff_abunds) != 0) {
        # remove seq from group totals
        if (private$has_group_data) {
          groups <- self$get_groups()
          for (j in seq_along(groups)) {
            private$group_totals[[groups[j]]] <-
              private$group_totals[[groups[j]]] - diff_abunds[j]
          }
        }
      }
    }
  )
)
