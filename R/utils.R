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
#' @import stringi
#'
#' @examples
#' (split_white_space("This is my string to split."))
#'
#' @returns A vector of Strings
split_white_space <- function(line) {
  words <- stringi::stri_split_charclass(line, "\\p{WHITE_SPACE}",
    omit_empty = TRUE
  )
  return(unlist(words))
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
