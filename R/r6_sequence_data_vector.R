#' @title sequence_data_vector
#' @description 'sequence_data_vector' is an R6 class that represents sequence
#'  FASTA and abundance data. It inherits from 'sequence_dataset'.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @import cli
#' @import r2r
#' @import data.table
#' @export
sequence_data_vector <- R6Class("sequence_data_vector",
  inherit = sequence_dataset,
  public = list(
    #' @description
    #' Create a new FASTA sequence dataset
    #' @param filename Name of file containing FASTA reads
    #' @param path Path to file containing FASTA reads
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
    #'     rchime_example("test.fasta"))
    #'
    #' @return A new `sequence_data_vector` object.
    initialize = function(filename = NULL, path = NULL) {
      super$initialize()
      private$count_table <- sequence_abundance_data$new()
      self$clear()
      if (!is.null(filename)) {
        df <- super$read_fasta_file(filename, path)
        private$names <- df$names
        private$sequences <- df$sequences
        private$comments <- df$comments
        private$trash_codes <- df$trash_codes
        private$starts <- df$starts
        private$ends <- df$ends
        private$lengths <- df$lengths
        private$ambigs <- df$ambigs
        private$polymers <- df$polymers
        private$numns <- df$numns

        num_table_rows <- 1
        all_names <- private$names
        for (i in seq_along(all_names)) {
          private$seq2row[all_names[i]] <- num_table_rows
          num_table_rows <- num_table_rows + 1
        }

        # all seqs are assummed 'good'
        rows <- length(private$names)
        private$table_seqs <- rep(TRUE, rows)

        abunds <- rep(1, rows)

        # add names with abundance 1
        private$count_table$add_seqs(names = private$names)
      }
      invisible(self)
    },

    #' @description
    #' Get summary of sequences data - (summary_seqs output)
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$print()
    print = function() {
      self$print_bad_accnos()
      summary <- summary_seqs(dataset = self, silent = TRUE)
      print(summary)
      private$count_table$print()
      cat(paste("\nNumber of unique seqs:", self$get_num_unique_seqs()), "\n")
      cat(paste("Total number of seqs:", self$get_num_seqs(), "\n"), "\n")
    },

    #' @description View summary of removed sequences
    #' 'print_bad_accnos()' outputs a summary of the reasons sequences were
    #'  removed from the dataset
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
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
    #'   dataset <- sequence_data_vector$new()
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
          (function(x) length(gregexpr("[N]", x)[[1]]))
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

        # # add new rows to dataset and set new seqs status to "good"
        private$names <- c(private$names, names)
        private$sequences <- c(private$sequences, sequences)
        private$comments <- c(private$comments, comments)
        private$trash_codes <- c(private$trash_codes, trash_codes)
        private$starts <- c(private$starts, starts)
        private$ends <- c(private$ends, ends)
        private$lengths <- c(private$lengths, lengths)
        private$ambigs <- c(private$ambigs, ambigs)
        private$polymers <- c(private$polymers, polymers)
        private$numns <- c(private$numns, numns)
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
            private$olengths <- c(private$olengths, olengths)
            private$ostarts <- c(private$ostarts, ostarts)
            private$oends <- c(private$oends, oends)
            private$mismatches <- c(private$mismatches, mismatches)
            private$ee <- c(private$ee, ee)
          } else {
            cli::cli_alert("The contigs data length is mismatched, removing.")
          }
        }

        # add names with abundance 1
        private$count_table$add_seqs(names = names)
      }

      invisible(self)
    },

    #' @description
    #' Remove all sequences from dataset
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$clear()
    clear = function() {
      # fasta data
      private$names <- c()
      private$sequences <- c()
      private$comments <- c()
      private$trash_codes <- c()
      private$starts <- c()
      private$ends <- c()
      private$lengths <- c()
      private$ambigs <- c()
      private$polymers <- c()
      private$numns <- c()

      # contigs data
      private$olengths <- c()
      private$ostarts <- c()
      private$oends <- c()
      private$mismatches <- c()
      private$ee <- c()
      # private$contigs_data <- data.table()

      # align data
      private$search_score <- c()
      private$sim_score <- c()
      private$longest_insert <- c()

      # abundance data
      private$count_table <- sequence_abundance_data$new()
      # map to quick reference sequence data
      private$seq2row <- hashmap()
      # map of reasons for removal
      private$accnos <- hashmap()

      # sequence counts
      private$total_bad <- 0
      private$unique_bad <- 0
      # private$total <- 0

      # included sequences
      private$table_seqs <- c()

      invisible(self)
    },

    #' @description
    #' Get abundance for sequence in the dataset
    #' @param name String, Name of sequence
    #' @param group (optional) Name of group to get abundance for,
    #'   if not provided the total abundance for the sequence is returned
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
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
        return(private$count_table$get_abund(name, group))
      } else {
        super$alert_eliminated_name(name, 0)
      }
      return(0)
    },

    #' @description
    #' Get abundances for sequence by group
    #' @param name String, Name of sequence
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   dataset$get_abunds(name =
    #'    "M00967_43_000000000-A3JHG_1_1111_8697_7063")
    #' @return Vector of abundances, in the same order as the groups
    #'             returned by dataset$get_groups().
    get_abunds = function(name) {
      return(private$count_table$get_abunds(name))
    },

    #' @description
    #' Get data.table containing a report of sequences removed
    #' @examples
    #'  dataset <- sequence_data_vector$new(filename =
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
    #'   dataset <- sequence_data_vector$new(filename =
    #'     rchime_example("test.fasta"))
    #'   reference <- sequence_data_vector$new(filename =
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
        searches <- private$search_score[private$table_seqs]
        sims <- private$sim_score[private$table_seqs]
        inserts <- private$longest_insert[private$table_seqs]

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
        lengths <- private$lengths[private$table_seqs]
        olengths <- private$olengths[private$table_seqs]
        ostarts <- private$ostarts[private$table_seqs]
        oends <- private$oends[private$table_seqs]
        mismatches <- private$mismatches[private$table_seqs]
        numns <- private$numns[private$table_seqs]
        ee <- private$ee[private$table_seqs]

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
    #'   dataset <- sequence_data_vector$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   count_table <- dataset$get_count_table()
    #' @return a data.table containing the count table data
    get_count_table = function() {
      df <- data.table(names = self$get_names())
      groups <- c()
      if (private$count_table$has_groups()) {
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
        df <- cbind(df, self$get_seqs_abunds(self$get_names()))
        names(df) <- c("names", "total")
      }
      return(df)
    },

    #' @description
    #' Get names of groups in the dataset
    #' @param name The name of the sequence you want groups for, optional
    #' @return A character vector
    get_groups = function(name = NULL) {
      return(private$count_table$get_groups(name))
    },

    #' @description
    #' Get the number of sequences represented in samples
    #' @param group String, name of sample
    #' @return A vector of integers
    get_group_totals = function(group = NULL) {
      return(private$count_table$get_group_totals(group))
    },

    #' @description
    #' Get names of sequences in the dataset
    #' @param group String, name of sample
    #' @examples
    #'   dataset <- sequence_data_vector$new(
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
      names <- private$names[private$table_seqs]

      if (!is.null(group)) {
        if (private$count_table$has_group(group)) {
          names_in_group <- c()
          for (name in names) {
            if (private$count_table$has_group(group, name)) {
              names_in_group <- c(names_in_group, name)
            }
          }
          return(names_in_group)
        }
      }
      return(names)
    },

    #' @description
    #' Get number of groups in the dataset
    #' @return An integer
    get_num_groups = function() {
      return(private$count_table$get_num_groups())
    },

    #' @description
    #' Get the number of sequences in the dataset
    #' @param group The name of the group you want number of sequences for,
    #'  optional
    #' @return An integer
    get_num_seqs = function(group = NULL) {
      if (is.null(group)) {
        return(private$count_table$get_total())
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
    #' Get data.table containing dataset
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   data <- dataset$export()
    #' @return data.table
    export = function() {
      data <- data.table(
        names = private$names,
        sequences = private$sequences,
        lengths = private$lengths,
        starts = private$starts,
        ends = private$ends,
        ambigs = private$ambigs,
        polymers = private$polymers,
        numns = private$numns
      )

      if (self$has_contigs_report()) {
        contigs_report <- data.table(
          olengths = private$olengths,
          ostarts = private$ostarts,
          oends = private$oends,
          mismatches = private$mismatches,
          ee = private$ee
        )

        data <- cbind(data, contigs_report)
      }

      if (self$has_align_report()) {
        align_report <- data.table(
          search_scores = private$search_score,
          sim_scores = private$sim_score,
          longest_inserts = private$longest_insert
        )

        data <- cbind(data, align_report)
      }

      data <- cbind(data, private$count_table$export())

      data <- cbind(data, data.table(
        comments = private$comments,
        trash_codes = private$trash_codes
      ))

      return(data)
    },

    #' @description
    #' Get vector containing FASTA nucleotide strings
    #' @param group String, name of sample
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
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
      if (!is.null(group)) {
        if (private$count_table$has_group(group)) {
          # names and sequences in table
          names <- self$get_names()
          seqs <- private$sequences[private$table_seqs]

          seqs_in_group <- c()
          for (i in seq_along(names)) {
            if (private$count_table$has_group(group, names[i])) {
              seqs_in_group <- c(seqs_in_group, seqs[i])
            }
          }
          return(seqs_in_group)
        }
      }

      return(private$sequences[private$table_seqs])
    },

    #' @description
    #' Get abundances for sequences in dataset
    #' @param names vector of Strings for which you want abundances for
    #' @param bysample Boolean, if true then abundance broken down by sample,
    #'                 if false, then total abundance is returned
    #' @examples
    #'   dataset <- sequence_data_vector$new(
    #'   filename = rchime_example('test.fasta'))
    #'
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'
    #'   # get abundances of all sequences in your dataset
    #'   abunds <- dataset$get_seqs_abunds()
    #'
    #'   # get abundances of sequences divided by sample
    #'   abunds <- dataset$get_seqs_abunds(bysample = TRUE)
    #' @return Vector of abundances, in the same order as the sequences
    #'             returned by dataset$get_names().
    get_seqs_abunds = function(bysample = FALSE) {
      names <- self$get_names()
      return(private$count_table$get_seqs_abunds(names, bysample))
    },

    #' @description
    #' Get data.table containing summary report data
    #' @examples
    #'   dataset <- sequence_data_vector$new(filename =
    #'   rchime_example("test.fasta"))
    #'   dataset <- dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'   df <- dataset$get_summary_report()
    get_summary_report = function() {
      lengths <- private$lengths[private$table_seqs]
      ambigs <- private$ambigs[private$table_seqs]
      polymers <- private$polymers[private$table_seqs]
      numns <- private$numns[private$table_seqs]
      starts <- private$starts[private$table_seqs]
      ends <- private$ends[private$table_seqs]

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
      if (length(private$search_score) == 0) {
        return(FALSE)
      }
      return(TRUE)
    },

    #' @description
    #' Determine if the dataset contains contigs report data
    #' @return Boolean
    has_contigs_report = function() {
      if (length(private$olengths) == 0) {
        return(FALSE)
      }
      return(TRUE)
    },

    #' @description
    #' Determine if a group is present in the dataset
    #' @param group String, Name of sample
    #' @return Boolean
    has_group = function(group) {
      if (private$count_table$has_group(group)) {
        return(TRUE)
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
    #'   dataset <- sequence_data_vector$new(filename =
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
    #'   dataset <- sequence_data_vector$new(filename =
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
          seqs_with_trashcode <- grepl(
            trash_code,
            private$trash_codes
          )

          # not all seqs with this trash code will be reinstated. The
          # sequence may have more than one trashcode
          names_2add <- private$get_seqs_for_reinstate(trash_code)

          accnos_names <- private$names[seqs_with_trashcode]

          for (j in seq_along(names_2add)) {
            seq_row <- private$seq2row[[names_2add[j]]]

            # set include seq flags
            private$table_seqs[seq_row] <- TRUE
            private$trash_codes[seq_row] <- ""

            # update groups and abunds in count_table
            private$count_table$reinstate_seq(names_2add[j])

            abund <- self$get_abund(names_2add[j])

            # update totals
            # private$total <- private$total + abund
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
    #'   dataset <- sequence_data_vector$new(filename =
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
          for (i in seq_along(names)) {
            row_to_remove <- private$seq2row[[names[i]]]
            if (!is.null(row_to_remove)) {
              private$remove_seq(names[i], trash_codes[i], TRUE)
            }
          }

          # private$total <- private$count_table$get_total()
        }
      }
      invisible(self)
    },

    #' @description
    #' Add group assignment data
    #' @param names a vector of sequence names
    #' @param groups a vector of group assignments
    #' @param abundances a vector of sample abundances
    #' @param filename String, Name of mothur formatted count file
    #' @param path String, Path to mothur formatted count file
    #' @examples
    #'
    #'  # mothur count file
    #'  # Representative_Sequence     total   sample2	sample3	sample4
    #'  # seq1	1150	250	400	500
    #'  # seq2	115	25	40	50
    #'  # seq3	50	25	25	0
    #'  # seq4	4	0	0	4
    #'
    #' # inputted as a sample table
    #' names <- c("seq1", "seq1", "seq1",
    #'           "seq2", "seq2", "seq2",
    #'           "seq3", "seq3",
    #'           "seq4")
    #' groups <- c("sample2", "sample3", "sample4",
    #'            "sample2", "sample3", "sample4",
    #'            "sample2", "sample3",
    #'            "sample4")
    #' abundances <- c(250, 400, 500,
    #'                25, 40, 50,
    #'                25, 25,
    #'                4)
    #'
    #' dataset <- sequence_data_vector$new()
    #' unique_names <- unique(names)
    #' sequences <- c("ATGGGCT", "..TG--ACCGT..", "..GGuatgc..", "..GGTAC-T..")
    #' dataset$add_seqs(unique_names, sequences)
    #' dataset$set_group_assignments(names, groups, abundances)
    #'
    #' # or read fasta file and mothur formatted count file
    #'
    #' dataset <- sequence_data_vector$new(filename =
    #'         rchime_example("test.fasta"))
    #' dataset$set_group_assignments(filename =
    #'         rchime_example("test.count_table"))
    #'
    set_group_assignments = function(names = NULL, groups = NULL,
                                     abundances = NULL,
                                     filename = NULL, path = NULL) {
      if (!is.null(names) && !is.null(groups)) {
        # find unique names
        unique_names <- unique(names)

        if (length(unique_names) != length(private$names)) {
          super$abort_length_mismatch(
            "names", "groups", length(unique_names),
            length(private$names)
          )
        }

        for (i in seq_along(unique_names)) {
          seq_info <- self$in_table(unique_names[i])
          if (!seq_info$include) {
            super$abort_name_missing(unique_names[i])
          }
        }

        private$count_table$set_group_assignments(names, groups, abundances)
      } else if (!is.null(filename)) {
        private$count_table$set_group_assignments(
          self$get_names(), NULL, NULL,
          filename, path
        )
      }

      invisible(self)
    },

    #' @description
    #' Set abundance data for sequences
    #' @param names a vector of sequence names
    #' @param abunds a list of abundance assignments
    #' @param reason, String (optional) reason for setting sequences to
    #'  zero abundance
    #' @examples
    #'   dataset <- sequence_data_vector$new()
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
              private$count_table$set_abundance(names[i], abunds[i])
            }
          }
        }
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
    #'   dataset <- sequence_data_vector$new()
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
          (function(x) length(gregexpr("[N]", x)[[1]]))
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

          private$sequences[indexes[i]] <- sequences[i]
          private$starts[indexes[i]] <- starts[i]
          private$ends[indexes[i]] <- ends[i]
          private$lengths[indexes[i]] <- lengths[i]
          private$ambigs[indexes[i]] <- ambigs[i]
          private$polymers[indexes[i]] <- polymers[i]
          private$numns[indexes[i]] <- numns[i]

          if (length(comments) != 0) {
            private$comments[indexes[i]] <- comments[i]
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

              private$search_score <- rep(0, num_all_seqs)
              private$sim_score <- rep(0, num_all_seqs)
              private$longest_insert <- rep(0, num_all_seqs)
            }

            # fill in "good" seqs align data
            for (i in seq_along(indexes)) {
              cli_progress_update()

              private$search_score[indexes[i]] <- search_scores[i]
              private$sim_score[indexes[i]] <- sim_scores[i]
              private$longest_insert[indexes[i]] <- longest_inserts[i]
            }
            cli_progress_done()
          } else {
            cli::cli_alert("The align data length is mismatched, removing.")
          }
        }
      }
      invisible(self)
    },

    #' @description
    #' Write count table data to mothur formatted count file
    #' @param filename String, name of file you want to write to
    #' @param path path to count file (optional)
    write_count_file = function(filename, path = NULL) {
      private$count_table$write_count_file(filename, path, self$get_names())
    }
  ),
  private = list(

    # FASTA data
    # names = c(), sequences = c(), comments = c(),
    # trash_codes = c(), starts = c(), ends = c(),
    # lengths = c(), ambigs = c(), polymers = c(), numns = c(),
    names = c(), sequences = c(), comments = c(),
    trash_codes = c(), starts = c(), ends = c(),
    lengths = c(), ambigs = c(), polymers = c(), numns = c(),
    count_table = NULL,

    # contigs report data
    # olengths = c(), ostarts = c(), oends = c(),
    # mismatches = c(), ee = c() same order as seq_data
    olengths = c(), ostarts = c(), oends = c(),
    mismatches = c(), ee = c(),

    # align report data
    # search_score, sim_score and longest_insert
    search_score = c(), sim_score = c(), longest_insert = c(),

    # sequence name -> row in seq_data
    # seq2row["seq1"] will return row in seq_data where seq1's info is
    seq2row = hashmap(),

    # vector of booleans indicating whether sequence is present in 'current'
    # table.
    table_seqs = c(),

    # trash_code -> count of seqs
    accnos = hashmap(),
    total_bad = 0,
    unique_bad = 0,
    total = 0,

    # Clear sequences from dataset
    finalize = function() {
      self$clear()
    },
    get_seqs_lengths = function() {
      return(unlist(lapply(
        private$sequences,
        (function(x) nchar(x))
      )))
    },

    # returns vector of names to add back into active state
    get_seqs_for_reinstate = function(trash_code) {
      # remove trash code from that sequences list of trash codes
      private$trash_codes <-
        gsub(
          trash_code, "",
          private$trash_codes
        )

      # remove double comma
      private$trash_codes <- gsub(",,", ",", private$trash_codes)

      # if it's now good, add the seq back in
      seqs_to_re_add <- (private$trash_codes == ",")

      # get the names
      reinstated_names <- private$names[seqs_to_re_add]
    },
    merge_seq = function(names, reason, group = NULL) {
      # make sure seqs are in the dataset
      for (j in seq_along(names)) {
        seq_index <- private$seq2row[[names[j]]]
        if (is.null(seq_index)) {
          super$abort_name_missing(names[j])
        }
      }

      if (length(names) != 1) {
        reason <- paste(reason, ",", sep = "")

        private$count_table$merge_seqs(names, group)
        names <- names[-c(1)]

        if (is.null(group)) {
          # remove duplicate names from "in table" status
          for (i in seq_along(names)) {
            private$remove_seq(names[i], reason, FALSE)
          }
        } else {
          for (i in seq_along(names)) {
            if (self$get_abund(names[i]) == 0) {
              private$remove_seq(names[i], reason, FALSE)
            }
          }
        }
      }

      invisible(self)
    },
    remove_seq = function(name, trash_code, update_totals = TRUE) {
      # when merging sequences we don't want to update group totals
      abunds <- self$get_abunds(name)
      abund <- sum(abunds)

      # if we are merging seqs, don't update group totals, or numSeqs
      if (update_totals) {
        private$count_table$remove_seq(name)
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
      private$trash_codes[row_to_remove] <- trash_code
      private$table_seqs[row_to_remove] <- FALSE

      invisible(self)
    }
  )
)
