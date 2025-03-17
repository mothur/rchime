#' @title sequence
#' @description 'sequence' is an R6 class representing a single FASTA
#' sequence read
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom waldo compare
#' @import cli
#' @import stringi
#' @export
sequence <- R6Class("sequence",
  public = list(

    #' @description
    #' Create a new sequence object.
    #' @param name Name of sequence
    #' @param sequence Sequence bases
    #' @param comment Sequence comment, optional
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' @return A new `sequence` object.
    initialize = function(name = "", sequence = "", comment = "") {
      self$set_name(name)
      self$set_sequence(sequence)
      self$set_comment(comment)
    },
    #' @description Print a FASTA sequence
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$print()
    print = function() {
      print_name <- paste(">", private$name, sep = "")
      cat(paste(print_name, private$comment[1], sep = "\t"),
        private$aligned,
        sep = "\n"
      )
    },

    #' @description
    #' Get string containing sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$get_inLine_seq()
    get_inLine_seq = function() {
      print_name <- paste(">", private$name, sep = "")
      output <- paste(paste(print_name, private$comment[1], sep = "\t"),
        private$aligned,
        sep = "\n"
      )
      return(output)
    },

    #' @description
    #' Get FASTA sequence name
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$get_name()
    get_name = function() {
      return(private$name)
    },

    #' @description
    #' Set FASTA sequence name
    #' @param name New sequence name.
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$set_name("newName")
    set_name = function(name) {
      private$name <- name
      invisible(self)
    },

    #' @description
    #' Set FASTA sequence data
    #' @param sequence New sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$set_sequence("--atgct--tt--ggcat")
    set_sequence = function(sequence = "") {
      # '--atgcu--tu--' becomes '..ATGCT--T..'

      if (sequence == "") {
        private$aligned <- ""
      } else {
        # RNA -> DNA and force uppercase
        private$aligned <- stri_replace_all_fixed(
          stri_trans_toupper(sequence), "U", "T"
        )

        # update gap char if needed
        length <- stri_length(private$aligned)
        first_non_gap <- stri_locate_first_regex(private$aligned, "[*.A-Z]")
        last_non_gap <- stri_locate_last_regex(private$aligned, "[*.A-Z]")

        # leading gaps
        if (first_non_gap[1] != 1) {
          replacement <- paste(rep(".", first_non_gap[1] - 1), collapse = "")
          private$aligned <- stri_sub_replace(private$aligned, 1,
            first_non_gap[1] - 1,
            value = replacement
          )
        }

        # trailing gaps
        if (last_non_gap[1] != length) {
          replacement <- paste(rep(".", (length - last_non_gap[1])),
            collapse = ""
          )
          private$aligned <-
            stri_sub_replace(private$aligned, last_non_gap[1] + 1, length,
              value = replacement
            )
        }
      }
      invisible(self)
    },

    #' @description
    #' Get aligned FASTA sequence data, if provided
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$get_aligned()
    get_aligned = function() {
      return(private$aligned)
    },

    #' @description
    #' Get unaligned FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$get_unaligned()
    get_unaligned = function() {
      return(stri_replace_all_regex(private$aligned, "[.-]", ""))
    },

    #' @description
    #' Set FASTA sequence comment
    #' @param comment New sequence comment
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$set_comment("This is my comment")
    set_comment = function(comment = "") {
      private$comment <- comment
      invisible(self)
    },

    #' @description
    #' Get FASTA sequence comment, if provided
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--tt--")
    #' fasta_read$get_comment()
    get_comment = function() {
      return(private$comment)
    },

    #' @description
    #' Get number of ambiguous bases in the FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_num_ambig()
    get_num_ambig = function() {
      length <- stri_length(private$aligned)
      num_non_ambig <- stri_count(private$aligned, regex = "[*.*-ATGC]")
      return(length - num_non_ambig)
    },

    #' @description
    #' Get length of the longest homopolymer in the FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_long_homop()
    get_long_homop = function() {
      longest_homop <- 1
      size <- 2

      # create homop vector of length 2 - start
      homops <- private$get_homop(size)

      # while you can still find homopolymers of this size
      while (any(stri_detect_fixed(private$aligned, homops))) {
        longest_homop <- size
        size <- size + 1

        # increase homop length
        homops <- private$get_homop(size)
      }
      return(longest_homop)
    },

    #' @description
    #' Get number of bases in the FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_num_bases()
    get_num_bases = function() {
      return(stri_length(stri_replace_all_regex(
        private$aligned,
        "[.-]", ""
      )))
    },

    #' @description
    #' Get the alignment position of the first base in the
    #'   FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_start()
    get_start = function() {
      return(stri_locate_first_regex(private$aligned, "[A-Z]")[1])
    },

    #' @description
    #' Get the alignment position of the last base in the
    #'   FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_end()
    get_end = function() {
      return(stri_locate_last_regex(private$aligned, "[A-Z]")[1])
    },

    #' @description
    #' Get the number of N's in the FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_num_ns()
    get_num_ns = function() {
      return(stri_count(private$aligned, regex = "[N]"))
    },

    #' @description
    #' Get the length of the alignment of the FASTA sequence data
    #' @examples
    #' fasta_read <- sequence$new("myseq", "--atgct--nn--")
    #' fasta_read$get_aligned_length()
    get_aligned_length = function() {
      return(stri_length(private$aligned))
    }
  ),
  private = list(
    aligned = "",
    comment = "",
    name = "",
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
