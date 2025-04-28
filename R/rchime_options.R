#' @title rchime_options
#'
#' @description
#' The rchime_options function allows you to set parameters for the
#' [[rchime()]] function. We recommend using the defaults, unless you
#' have a compelling reason to do so.
#'
#' "UCHIME parameters are optimized for detection of very low-divergence
#' chimeras. In typical applications such as 16S OTU picking from next-gen
#' reads, chimeras over divergence less than the OTU radius may not be
#' important, in which case it may be better to retune parameters.
#' This can be done by increasing ‑mindiv and reducing -minh."
#'
#' @references
#'  Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
#'  UCHIME improves sensitivity and speed of chimera detection.
#'  Bioinformatics 27:2194.
#'
#' @param abskew Float, the minimum abundance skew (denovo only). Default = 2.0.
#'               abskew <- min (abund(parent1), abund(parent2)) / abund(query)
#' @param minh Float, Mininum score to report chimera. Default 0.3. Values from
#'     0.1 to 5 might be reasonable. Lower values increase sensitivity
#'    but may report more false positives. If you decrease --xn,
#'    you may need to increase --minh, and vice versa.
#' @param mindiv Float, Minimum divergence ratio, default 0.5. Div ratio is
#'    100%% %%identity between query sequence and the closest candidate for
#'    being a parent. If you don't care about very close chimeras,
#'    then you could increase --mindiv to, say, 1.0 or 2.0, and
#'    also decrease --min h, say to 0.1, to increase sensitivity.
#'    How well this works will depend on your data. Best is to
#'    tune parameters on a good benchmark.
#' @param xn Float, Weight of a no vote, also called the beta parameter.
#'    Default 8.0. Decreasing this weight to around 3 or 4 may give better
#'    performance on denoised data.
#' @param dn Float, Pseudo-count prior on number of no votes. Default 1.4.
#'    Probably no good reason to change this unless you can retune to a good
#'    benchmark for your data. Reasonable values are probably in the
#'    range from 0.2 to 2.
#' @param xa Float, Weight of an abstain vote. Default 1. So far, results do not
#'    seem to be very sensitive to this parameter, but if you have
#'    a good training set might be worth trying. Reasonable values
#'    might range from 0.1 to 2.
#' @param chunks Integer, Number of chunks to extract from the query sequence
#'    when searching for parents. Default 4.
#' @param minchunk Integer, Minimum length of a chunk. Default 64.
#' @param idsmoothwindow Integer, Length of id smoothing window. Default 32.
#' @param maxp Integer, Maximum number of candidate parents to consider.
#'   Default 2. In tests so far, increasing maxp gives only a very small
#'   improvement in sensivity but tends to increase the error rate quite a bit.
#' @param skipgaps, bool, If skipgaps is specified, columns containing gaps do
#'  not count as diffs. Default = TRUE.
#' @param skipgaps2 bool, If skipgaps2 is specified, if column is immediately
#'  adjacent to a column containing a gap, it is not counted as a diff.
#'  Default = TRUE.
#' @param minlen Integer, Minimum sequence length. Default = 10.
#' @param maxlen Integer, Maximum sequence length. Default = 10000.
#'
#' @return List containing selected rchime parameters and their values.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @export
#'
rchime_options <- function(
    abskew = 2.0,
    minh = 0.3, mindiv = 0.5, xn = 8.0,
    dn = 1.4, xa = 1.0, chunks = 4, minchunk = 64,
    idsmoothwindow = 32, maxp = 2, skipgaps = TRUE,
    skipgaps2 = TRUE, minlen = 10, maxlen = 10000) {
  parameters <- list()

  # if the user sets uchime options then add to parameters
  if (abskew != 2.0) {
    parameters <- c(parameters, abskew = abskew)
  }
  if (minh != 0.3) {
    parameters <- c(parameters, minh = minh)
  }
  if (mindiv != 0.5) {
    parameters <- c(parameters, mindiv = mindiv)
  }
  if (xn != 8.0) {
    parameters <- c(parameters, xn = xn)
  }
  if (dn != 1.4) {
    parameters <- c(parameters, dn = dn)
  }
  if (xa != 1.0) {
    parameters <- c(parameters, xa = xa)
  }
  if (chunks != 4) {
    parameters <- c(parameters, chunks = chunks)
  }
  if (minchunk != 64) {
    parameters <- c(parameters, minchunk = minchunk)
  }
  if (idsmoothwindow != 32) {
    parameters <- c(parameters, idsmoothwindow = idsmoothwindow)
  }
  if (maxp != 2) {
    parameters <- c(parameters, maxp = maxp)
  }
  if (!skipgaps) {
    parameters <- c(parameters, skipgaps = FALSE)
  }
  if (!skipgaps2) {
    parameters <- c(parameters, skipgaps2 = FALSE)
  }
  if (minlen != 10) {
    parameters <- c(parameters, minlen = minlen)
  }
  if (maxlen != 10000) {
    parameters <- c(parameters, maxlen = maxlen)
  }

  return(parameters)
}
