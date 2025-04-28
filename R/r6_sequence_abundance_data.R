#' @title sequence_abundance_data
#' @description 'sequence_abundance_data' is an R6 class that stores the
#'  abundance data of sequences by sample.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @import cli
#' @import r2r
#' @export
sequence_abundance_data <- R6Class("sequence_abundance_data",
  inherit = sequence_dataset,
  public = list(
    #' @description
    #' Create a new abundance table
    #' @return A new `sequence_abundance_data` object.
    initialize = function() {
      super$initialize()
      self$clear()
      invisible(self)
    },

    #' @description View summary of count table data
    #' 'print()' outputs a summary of the sequence abundance data
    #' @examples
    #'   counts <- sequence_abundance_data$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   groups <- c("sample1", "sample2", "sample1")
    #'   counts$set_group_assignments(names, groups)
    #'   counts$print()
    print = function() {
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
    #' Adds sequences with abundance of 1
    #' @param names a vector of sequence names
    #' @examples
    #'
    #'   dataset <- sequence_abundance_data$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   dataset$add_seqs(names)
    #'
    add_seqs = function(names) {
      if (private$has_group_data) {
        cli::cli_abort("[ERROR]: The 'sequence_data_vector' contains group
                      information, you must provide by sample data.")
      } else {
        # provide total abundance for each seq rather than by sample
        for (i in seq_along(names)) {
          private$counts[[names[i]]] <- list(1)
        }
        private$table_groups <- c(TRUE)
        private$total <- private$total + length(names)
      }
      invisible(self)
    },

    #' @description
    #' Remove all sequences from dataset
    #' @examples
    #'   dataset <- sequence_abundance_data$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   dataset$add_seqs(names)
    #'   dataset
    #'   dataset$clear()
    #'   dataset
    clear = function() {
      private$counts <- hashmap()

      # group names mapped to column in counts
      private$groups <- hashmap()

      # map of pre-calculated group totals
      private$group_totals <- hashmap()

      private$total <- 0
      private$has_group_data <- FALSE

      invisible(self)
    },

    #' @description
    #' Get data.table containing abundance data
    #' @param names vector names to export (optional)
    #' @return data.table
    export = function(names = NULL) {
      # use all names in counts
      if (is.null(names)) {
        names <- unlist(keys(private$counts))
      }

      abunds <- self$get_seqs_abunds(names, private$has_group_data)

      if (private$has_group_data) {
        groups <- self$get_groups()

        total_abundance <- rep(0, length(abunds))
        for (i in seq_along(abunds)) {
          total_abundance[i] <- sum(abunds[[i]])
        }

        data <- data.table(total_abundance = total_abundance)
        for (group in groups) {
          group_index <- private$groups[[group]]

          group_abunds <- rep(0, length(abunds))
          for (i in seq_along(abunds)) {
            group_abunds[i] <- abunds[[i]][group_index]
          }
          data <- cbind(data, group_abunds)
        }

        names(data) <- c("total_abundance", groups)
        return(data)
      } else {
        data <- data.table(total_abundance = abunds)
        return(data)
      }
      return(data.table())
    },

    #' @description
    #' Get abundance for sequence in the dataset
    #' @param name String, Name of sequence
    #' @param group (optional) Name of group to get abundance for,
    #'   if not provided the total abundance for the sequence is returned
    #' @return An integer
    get_abund = function(name, group = NULL) {
      data <- private$counts[[name]]

      # if the name is NOT in the table
      if (is.null(data)) {
        cli::cli_abort("{.var{name}} is not present in your dataset, aborting.")
      }

      if (is.null(group)) {
        if (private$has_group_data) {
          return(sum(self$get_abunds(name)))
        } else {
          return(private$counts[[name]][[1]][1])
        }
      } else if (self$has_group(group)) {
        group_index <- private$groups[[group]]
        seqs_inc_index <- which(data[[1]] %in% group_index)
        if (length(seqs_inc_index) != 0) {
          return(data[[2]][seqs_inc_index])
        }
      }
      return(0)
    },

    #' @description
    #' Get abundances for sequence by group
    #' @param name String, Name of sequence
    #' @return Vector of abundances, in the same order as the groups
    #'             returned by counts$get_groups().
    get_abunds = function(name) {
      data <- private$counts[[name]]

      # if the name is in the table
      if (!is.null(data)) {
        if (private$has_group_data) {
          # c(3,2) c(100, 45)
          abunds <- rep(0, length(private$groups))
          # group_abund_data[[1]] -> c(3,2)
          for (i in seq_along(data[[1]])) {
            # group_abund_data[[2]] -> c(100,45)
            # i -> 1, abunds[3] <- 100
            # i -> 2, abunds[2] <- 45
            abunds[data[[1]][i]] <- data[[2]][i]
          }
          return(abunds[private$table_groups])
        } else {
          return(private$counts[[name]][[1]])
        }
      } else {
        # TODO ERROR
        return(c())
      }
    },

    #' @description
    #' Get names of groups in the dataset
    #' @param name The name of the sequence you want groups for, optional
    #' @return A character vector
    get_groups = function(name = NULL) {
      if (is.null(name)) {
        return(private$get_included_groups())
      } else {
        data <- private$counts[[name]]
        if (!is.null(data) && private$has_group_data) {
          inc_groups <- private$get_included_groups()
          included_indexes <- private$get_included_groups_indexes()
          seqs_group_inds <- which(data[[1]] %in% included_indexes)
          # returns groups that are included in dataset and present in seq
          return(inc_groups[data[[1]][seqs_group_inds]])
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
    #' Get number of groups in the dataset
    #' @return An integer
    get_num_groups = function() {
      if (private$has_group_data) {
        return(sum(private$table_groups, na.rm = TRUE))
      }
      return(0)
    },

    #' @description
    #' Get sample table returns data.table containing the 3 columns: name,
    #' group and abundance.
    #' @param names, vector of names to include in table (optional)
    #' @examples
    #'   dataset <- sequence_abundance_data$new()
    #'   dataset$set_group_assignments(filename =
    #'   rchime_example("test.count_table"))
    #'
    #'   sample_table <- dataset$get_sample_table()
    #'
    #' @return data.table
    get_sample_table = function(names = NULL) {
      # use all names in counts
      if (is.null(names)) {
        names <- unlist(keys(private$counts))
      }

      # if groups, create a row for each seq,sample,abund
      if (private$has_group_data) {
        table_names <- c()
        table_groups <- c()
        table_abundances <- c()

        groups <- self$get_groups()

        for (name in names) {
          data <- private$counts[[name]]

          table_names <- c(table_names, rep(name, length(data[[1]])))
          # add each group name is present in
          for (i in seq_along(data[[1]])) {
            table_groups <- c(table_groups, groups[data[[1]][i]])
          }
          # add each abundance for each sample name is present in
          for (i in seq_along(data[[2]])) {
            table_abundances <- c(table_abundances, data[[2]][i])
          }
        }
        data <- data.table(
          name = table_names, group = table_groups,
          abundance = table_abundances
        )
        return(data)
      } else {
        abunds <- self$get_seqs_abunds(names)
        data <- data.table(name = names, abundance = abunds)
        return(data)
      }
      return(data.table())
    },

    #' @description
    #' Get number of sequences in the dataset
    #' @return An integer
    get_total = function() {
      return(private$total)
    },

    #' @description
    #' Get abundances for sequences in dataset
    #' @param names vector of Strings for which you want abundances for
    #' @param bysample Boolean, if true then abundance broken down by sample,
    #'                 if false, then total abundance is returned
    #' @return Vector of abundances, in the same order as the sequences
    #'             in names.
    get_seqs_abunds = function(names, bysample = FALSE) {
      if (bysample) {
        abunds <- list(length(names))
        for (i in seq_along(names)) {
          abunds[[i]] <- self$get_abunds(names[i])
        }
        return(abunds)
      } else {
        abunds <- rep(0, length(names))
        for (i in seq_along(names)) {
          abunds[i] <- self$get_abund(names[i])
        }
        return(abunds)
      }
      return(c())
    },

    #' @description
    #' Determine if the abundance data is divided by sample
    #' @return Boolean
    has_groups = function() {
      return(private$has_group_data)
    },

    #' @description
    #' Determine if a group is present in the dataset
    #' @param group String, Name of sample
    #' @param name String, name of sequence. If name is provided, determine if
    #' sequence is present in the sample (optional)
    #' @return Boolean
    has_group = function(group, name = NULL) {
      if (private$has_group_data) {
        all_groups <- self$get_groups(name)
        return(group %in% all_groups)
      }
      return(FALSE)
    },

    #' @description merges identical sequences in the dataset
    #' @param names Vector of names to be merged
    #' @param group String, when a group is provided only that groups counts
    #'                       will be merged (optional)
    #' @param reason String, reason for merging sequence
    #' @examples
    #'   dataset <- sequence_abundance_data$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   dataset$add_seqs(names)
    #'   dataset
    #'   dataset$merge_seqs(c("seq1", "seq2"))
    #'   dataset
    merge_seqs = function(names, group = NULL) {
      if (length(names) != 1) {
        group_index <- private$groups[[group]]
        keeper_name <- names[1]
        names <- names[-c(1)]

        k_abunds <- self$get_abunds(keeper_name)

        for (i in seq_along(names)) {
          dup_name <- names[i]
          dup_abunds <- self$get_abunds(dup_name)

          if (is.null(group_index)) {
            # merge all counts
            k_abunds <- k_abunds + dup_abunds
          } else {
            # merge counts for this group
            k_abunds[group_index] <- k_abunds[group_index] +
              dup_abunds[group_index]
            dup_abunds[group_index] <- 0
            d <- list(dup_abunds)
            private$counts[[dup_name]] <- private$to_sparse(d)
          }
        }

        # add dups groups and abunds to keeper
        k_abunds <- list(k_abunds)
        private$counts[[keeper_name]] <- private$to_sparse(k_abunds)
      }

      invisible(self)
    },

    #' @description
    #' Reinstates sequence in active table
    #' @param name name of sequence to reinstate
    reinstate_seq = function(name) {
      groups <- private$get_included_groups()

      # check if this seq re-adds a group
      this_seqs_groups <- private$get_groups_for_reinstate(name)

      groups_to_re_add <- setdiff(this_seqs_groups, groups)

      if (length(groups_to_re_add) != 0) {
        for (k in seq_along(groups_to_re_add)) {
          group_index <- private$groups[[groups_to_re_add[k]]]
          private$table_groups[group_index] <- TRUE
        }
        groups <- self$get_groups()
      }

      abunds <- self$get_abunds(name)

      private$total <- private$total + sum(abunds)

      if (private$has_group_data) {
        diff_abunds <- abunds * -1L
        private$update_group_totals(diff_abunds)
      }

      invisible(self)
    },

    #' @description
    #' Remove sequences from dataset
    #' @param name String, name of sequence to remove
    remove_seq = function(name) {
      # when merging sequences we don't want to update group totals
      abunds <- self$get_abunds(name)
      abund <- sum(abunds)

      # remove seq from group totals
      if (private$has_group_data) {
        private$update_group_totals(abunds)
      }
      private$total <- private$total - abund

      # are their samples to remove now? This can happen if the seqs
      # removed eliminate a group
      new_totals <- self$get_group_totals()

      # are there groups in the dataset with no sequences
      if (0 %in% new_totals) {
        groups <- self$get_groups()
        for (i in seq_along(groups)) {
          if (new_totals[i] == 0) {
            private$remove_group(groups[i])
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Set abundance data for sequence
    #' @param name a vector of sequence names
    #' @param abundance a list of abundance assignments
    #' @examples
    #'   dataset <- sequence_abundance_data$new()
    #'   names <- c("seq1", "seq2", "seq3")
    #'   dataset$add_seqs(names)
    #'
    #'   # set total abundance. seq1 -> 100, seq2 -> 50, seq3 -> 5
    #'   abunds <- list(100, 50, 5)
    #'   dataset$set_abundance("seq1", abunds[1])
    #'   dataset$set_abundance("seq2", abunds[2])
    #'   dataset$set_abundance("seq2", abunds[3])
    #'
    #'   # or with group data
    #'   dataset <- sequence_abundance_data$new()
    #'   groups <- c("sample1", "sample1", "sample2")
    #'   dataset$set_group_assignments(names, groups)
    #'
    #'   # set the abundance of 'seq1' in 'sample1' to 150,
    #'   #                      'seq1' in 'sample2' to 0
    #'   dataset$set_abundance("seq1", list(c(150, 0)))
    #'
    #'   # set the abundance of 'seq3' in 'sample1' to 75,
    #'   #                      'seq3' in 'sample2' to 5
    #'   dataset$set_abundance("seq3", list(c(75, 5)))
    #'   dataset
    set_abundance = function(name, abundance) {
      # save old abunds if you have groups so we can update group_totals
      orig_abunds <- self$get_abunds(name)

      # update counts
      private$counts[[name]] <- private$to_sparse(abundance)

      diff_abunds <- orig_abunds - abundance[[1]]
      private$total <- private$total - sum(diff_abunds)

      if (private$has_group_data) {
        private$update_group_totals(diff_abunds)
      }

      invisible(self)
    },

    #' @description
    #' Add group assignment data, clears out old sample data
    #' @param names a vector of sequence names
    #' @param groups a vector of group assignments
    #' @param abundances a vector of sample abundances
    #' @param filename String, Name of mothur formatted count file
    #' @param path String, Path to mothur formatted count file
    #' @examples
    #'  # mothur count file
    #'  # Representative_Sequence     total   sample2	sample3	sample4
    #'  # seq1	1150	250	400	500
    #'  # seq2	115	25	40	50
    #'  # seq3	50	25	25	0
    #'  # seq4	4	0	0	4
    #'
    #' # as a sample table
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
    #' abunds <- sequence_abundance_data$new()
    #' abunds$set_group_assignments(names, groups, abundances)
    #'
    #' # or read mothur formatted count file
    #'
    #' abunds <- sequence_abundance_data$new()
    #' abunds$set_group_assignments(filename =
    #'         rchime_example("test.count_table"))
    #'
    set_group_assignments = function(names = NULL, groups = NULL,
                                     abundances = NULL,
                                     filename = NULL, path = NULL) {
      self$clear()

      if (!is.null(names) && !is.null(groups)) {
        if (length(names) != length(groups)) {
          super$abort_length_mismatch(
            "names", "groups", length(names),
            length(groups)
          )
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

        sums <- rep(0, num_groups)

        if (num_groups != 0) {
          unique_names <- unique(names)

          # assume each abundance is 1 if not provided
          if (is.null(abundances)) {
            abundances <- rep(1, length(groups))
          }

          # assigning each sequence to a sample - make_contigs
          if (length(unique_names) == length(groups)) {
            for (i in seq_along(groups)) {
              group_index <- private$groups[[groups[i]]]
              private$counts[[names[i]]] <- list(
                c(group_index),
                abundances[i]
              )
              sums[group_index] <- sums[group_index] + abundances[i]
            }
          } else {
            # assigning sequences to multiple groups

            # set all counts entries to long format of abundance
            long_abundances <- rep(0, num_groups)
            for (name in unique_names) {
              private$counts[[name]] <- long_abundances
            }

            for (i in seq_along(groups)) {
              group_index <- private$groups[[groups[i]]]
              private$counts[[names[i]]][group_index] <- abundances[i]
              sums[group_index] <- sums[group_index] + abundances[i]
            }

            for (name in unique_names) {
              private$counts[[name]] <- private$to_sparse(
                list(private$counts[[name]])
              )
            }
          }

          # set group totals
          for (i in seq_along(unique_groups)) {
            private$group_totals[[unique_groups[i]]] <- sums[i]
          }
          private$total <- sum(sums)
        }
      } else if (!is.null(filename)) {
        private$read_count_file(names, filename, path)
      }
      invisible(self)
    },

    #' @description
    #' Write count table data to mothur formatted count file
    #' @param filename String, name of file you want to write to
    #' @param path path to count file (optional)
    #' @param fasta_names vector of strings containing names to write
    write_count_file = function(filename, path = NULL, fasta_names) {
      if (file.exists(filename)) {
        file.remove(filename)
      }

      file_conn <- file(filename, open = "a")
      group_names <- self$get_groups()

      writeLines(private$get_header(group_names), file_conn)

      if (private$has_group_data) {
        for (i in seq_along(fasta_names)) {
          abunds <- private$get_full_abunds(fasta_names[i])
          output <- paste(c(fasta_names[i], sum(abunds), abunds),
            collapse = "\t"
          )

          writeLines(output, file_conn)
        }
      } else {
        for (i in seq_along(fasta_names)) {
          output <- paste(c(fasta_names[i], self$get_abund(fasta_names[i])),
            collapse = "\t"
          )
          writeLines(output, file_conn)
        }
      }

      close(file_conn)
    }
  ),
  private = list(

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

    # groupName -> total abundance for group
    # sample4 -> 250 means sample4 contains 250 sequences
    group_totals = hashmap(),

    # vector of booleans indicating whether sample is present in 'current'
    # table.
    table_groups = c(),

    # is abundance data parsed by sample
    has_group_data = FALSE,
    total = 0,
    add_groups = function(groups) {
      groups <- sort(groups)
      private$table_groups <- c(rep(TRUE, length(groups)))

      for (i in seq_along(groups)) {
        private$groups[[groups[i]]] <- i
        private$group_totals[[groups[i]]] <- 0
      }
      private$has_group_data <- TRUE

      invisible(self)
    },

    # abunds a vector of abundances in full format, convert to sparse format
    # seq1 -> c(0, 5, 0, 2, 0, 3) becomes
    # seq1 -> list(c(4, 2, 6), c(2, 5, 3))
    # 2 from 'sample4', 5 from 'sample2', 3 from 'sample6'
    to_sparse = function(abunds) {
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

    # abunds a vector of abundances in full format
    # seq1 -> list(c(4, 2, 6), c(2, 5, 3)) becomes
    # seq1 -> c(0, 5, 0, 2, 0, 3)
    # 2 from 'sample4', 5 from 'sample2', 3 from 'sample6'
    get_full_abunds = function(name) {
      data <- private$counts[[name]]

      abunds <- rep(0, self$get_num_groups())

      for (i in seq_along(data[[1]])) {
        abunds[data[[1]][i]] <- data[[2]][i]
      }
      return(abunds)
    },

    # Compressed Format: groupIndex,abundance. 1,6 means an ...
    # 	2,sample2	3,sample3	1,sample4
    # Representative_Sequence     total  sample2	sample3	sample4
    # returns data.frame contain count tables abundance data
    fill_from_compressed_format = function(file_data, fasta_names) {
      # extract group names
      # line 2 looks like: "#2,sample2	3,sample3	1,sample4"
      # remove first '#'
      pieces <- strsplit(file_data[2], "#")[[1]]
      file_data[2] <- pieces[nzchar(pieces)]
      # file_data[2] <- stri_sub(file_data[2], from = comment[1] + 1)
      words <- split_white_space(file_data[2])
      num_seqs <- length(file_data) - 3

      has_group_data <- TRUE
      groups <- c()

      for (i in seq_along(words)) {
        # parse group name
        file_index <- split_at_char(words[i], ",")

        # save group names
        groups <- c(groups, file_index[2])
      }

      private$add_groups(groups)
      sums <- rep(0, length(private$table_groups))

      # read compressed data lines
      for (i in 4:length(file_data)) {
        words <- split_white_space(file_data[i])

        name <- words[1]

        # sanity check
        if (!is.null(fasta_names)) {
          if (!(name %in% fasta_names)) {
            cli::cli_abort("[ERROR]: The sequence {.var {name}} is in
                    your count file and not in your dataset.")
          }
        }

        abunds <- rep(0, length(words) - 2)
        for (j in 3:length(words)) {
          data <- split_at_char(words[j], ",")
          abunds[as.integer(data[1])] <- as.integer(data[2])
        }
        group_indexes <- which(abunds != 0)
        private$counts[[name]] <-
          list(group_indexes, abunds[group_indexes])
        sums[group_indexes] <- sums[group_indexes] + abunds[group_indexes]
      }

      # include all groups
      groups <- self$get_groups()

      for (i in seq_along(groups)) {
        private$group_totals[[groups[i]]] <- sums[i]
      }
      private$total <- sum(sums)

      invisible(self)
    },

    # returns data.frame contain count tables abundance data
    fill_from_uncompressed_format = function(file_data, fasta_names) {
      # uncompressed format
      # Representative_Sequence  total  sample2	sample3	sample4
      words <- split_white_space(file_data[1])

      num_seqs <- length(file_data) - 1
      private$has_group_data <- TRUE

      # no group data in file
      if (length(words) == 2) {
        private$table_groups <- c(TRUE)
        private$has_group_data <- FALSE
      } else {
        groups <- c()
        for (i in 3:length(words)) {
          groups <- c(groups, words[i])
        }

        private$add_groups(groups)
      }

      sums <- rep(0, length(private$table_groups))

      # read uncompressed data
      for (i in 2:length(file_data)) {
        words <- split_white_space(file_data[i])

        if (length(words) >= 2) {
          name <- words[1]
          seq_total <- as.integer(words[2])

          # sanity check
          if (!is.null(fasta_names)) {
            if (!(name %in% fasta_names)) {
              cli::cli_abort("[ERROR]: The sequence {.var {name}} is in
                    your count file and not in your dataset.")
            }
          }

          if (private$has_group_data) {
            abunds <- rep(0, self$get_num_groups())
            for (j in 3:length(words)) {
              # add all samples
              abunds[j - 2] <- as.integer(words[j])
            }

            group_indexes <- which(abunds != 0)
            private$counts[[name]] <-
              list(group_indexes, abunds[group_indexes])
            sums[group_indexes] <- sums[group_indexes] + abunds[group_indexes]
          } else {
            private$counts[[name]] <- list(seq_total)
            private$total <- private$total + seq_total
          }
        }
      }

      # set group totals
      if (private$has_group_data) {
        # include all groups
        groups <- self$get_groups()

        for (i in seq_along(groups)) {
          private$group_totals[[groups[i]]] <- sums[i]
        }
        private$total <- sum(sums)
      }
      invisible(self)
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

    # Get names of groups in the dataset
    get_included_groups = function() {
      return(sort(unlist(keys(private$groups))[private$table_groups]))
    },
    get_included_groups_indexes = function() {
      return(which(private$table_groups))
    },
    get_groups_for_reinstate = function(name) {
      all_groups <- sort(unlist(keys(private$groups)))
      indexes <- private$counts[[name]][[1]]
      return(all_groups[indexes])
    },


    # Read a mothur formatted count file
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

      # check for group info in header
      if (num_lines > 1) {
        comment <- regexpr("#", file_data[2])

        # is this file compressed
        if (comment[1] != -1) {
          private$fill_from_compressed_format(file_data, names)
        } else {
          private$fill_from_uncompressed_format(file_data, names)
        }
      }
      invisible(self)
    },


    # Remove group from dataset
    # @param group String, name of sample to remove
    remove_group = function(group) {
      group_index <- private$groups[[group]]
      if (!is.null(group_index)) {
        private$table_groups[group_index] <- FALSE
      }
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
      invisible(self)
    },
    finalize = function() {}
  )
)
