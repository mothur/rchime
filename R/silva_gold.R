#' @title silva_gold
#' @name silva_gold
#' @rdname silva_gold
#'
#' @description
#' The `silva_gold()` function creates a silva based reference for use with the
#' `rchime()` function in reference mode. It contains 5,181 sequences.
#'
#' @examples
#'
#' reference <- silva_gold()
#'
#' @return data.frame containing
#'   \href{https://mothur.org/wiki/silva_reference_files/}{silva reference files
#'   for chimera detection}
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @export
silva_gold <- function() {
  readRDS(rchime_example("silva.gold.rds"))
}
