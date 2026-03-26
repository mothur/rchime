#' Get path to rchime example
#'
#' rchime comes bundled with some example files in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @export
#' @examples
#' rchime_example()
#' rchime_example("reference.rds")
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
