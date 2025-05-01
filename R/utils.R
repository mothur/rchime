# =========================================================================== #
#' @title remove_file
#' @description Remove file, if it exists
#' @param filename, String containing name of file to remove
remove_file <- function(filename) {
  if (file.exists(filename)) {
    file.remove(filename)
  }
}
# =========================================================================== #
#' @title split_white_space
#' @description Split string at white space
#' @param line, String containing data to split
#'
#' @examples
#' (split_white_space("This is my string to split."))
#'
#' @returns A vector of Strings
split_white_space <- function(line) {
  words <- strsplit(line, "\\s")[[1]]
  words <- words[nzchar(x = words)]
  return(words)
}
# =========================================================================== #
#' @title split_at_char
#' @description Split string by deliminating character
#' @param line, String containing data to split
#' @param delim, Character to split data by
#' @examples
#' (split_at_char("This is the first piece, and this is the second", ","))
#'
#' @returns A vector of Strings
# split at character
split_at_char <- function(line, delim) {
  words <- strsplit(line, delim)
  return(unlist(words))
}

# =========================================================================== #
#' @title extract_name
#' @description extracts name from fasta file entry
#' @param line, String containing name data
#' @examples
#' (extract_name("mySeq This is my comment"))
#'
#' @returns A string
extract_name = function(line) {
  # name. comment
  # mySeq This is my comment
  val <- regexpr("\\s", line)
  if (val != -1) {
    name_comment <- substring(line, c(1, val + 1), c(val - 1, nchar(line)))
    return(name_comment[1])
  }
  return(line)
}
# =========================================================================== #
