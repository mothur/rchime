#' Get path to rchime example files
#'
#' rchime comes bundled with some example files in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#'
#' @return string, full path to file requested
#' @export
#' @examples
#' rchime_example()
#' rchime_example("silva.gold.rds")
rchime_example <- function(file = NULL) {
  path <- ""

  if (is.null(file)) {
    path <- system.file("extdata", package = "rchime")
  } else {
    path <- system.file("extdata", file,
      package = "rchime",
      mustWork = TRUE
    )
  }
  return(path)
}

#' @keywords internal
fill_required_parameters <- function(df, given_column_name) {
  if (!(given_column_name %in% names(df))) {
    message <- paste0(
      "Expected a data.frame column named ",
      given_column_name, " to be provided."
    )

    cli::cli_abort(message)
  }

  return(df[[given_column_name]])
}
