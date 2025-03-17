#' Get path to mothur2 example
#'
#' mothur2 comes bundled with some example files in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @export
#' @examples
#' mothur2_example()
#' mothur2_example("test.fasta")
mothur2_example <- function(file = NULL) {
  path <- ""

  if (is.null(file)) {
    path <- system.file("extdata", package = "mothur2")
  } else {
    path <- system.file("extdata", file,
      package = "mothur2",
      mustWork = TRUE
    )
  }
  return(path)
}
