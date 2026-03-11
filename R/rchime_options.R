#' @title rchime_options
#'
#' @description
#' The rchime_options function allows you to set parameters for the
#' [[rchime()]] function. We recommend using the defaults, unless you
#' have a compelling reason to do otherwise.
#'
#' @references
#'   Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
#'   versatile open source tool for metagenomics. PeerJ 4:e2584.
#'   doi: 10.7717/peerj.2584
#'
#' @references
#'  Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
#'  UCHIME improves sensitivity and speed of chimera detection.
#'  Bioinformatics 27:2194.
#'
#' @param abskew Float, the minimum abundance skew (denovo only). Default = 2.0.
#'               abskew <- min (abund(parent1), abund(parent2)) / abund(query)
#' @param minh Float, Mininum score to report chimera. Default 0.28. Values from
#'     0.1 to 5 might be reasonable. Lower values increase sensitivity
#'    but may report more false positives. If you decrease --xn,
#'    you may need to increase --minh, and vice versa.
#' @param mindiv Float, Minimum divergence ratio, default 0.8. Div ratio is
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
#' @param maxp Integer, Maximum number of candidate parents to consider.
#'   Default 3.
#'
#' @return List containing selected rchime parameters and their values.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#' @export
#'
rchime_options <- function(
  abskew = 2.0,
  minh = 0.28, mindiv = 0.8, xn = 8.0,
  dn = 1.4, maxp = 3
) {
  parameters <- list()

  # if the user sets uchime options then add to parameters
  if (abskew != 2.0) {
    parameters <- c(parameters, abskew = abskew)
  }
  if (minh != 0.28) {
    parameters <- c(parameters, minh = minh)
  }
  if (mindiv != 0.8) {
    parameters <- c(parameters, mindiv = mindiv)
  }
  if (xn != 8.0) {
    parameters <- c(parameters, xn = xn)
  }
  if (dn != 1.4) {
    parameters <- c(parameters, dn = dn)
  }
  if (maxp != 3) {
    parameters <- c(parameters, maxp = maxp)
  }
  return(parameters)
}
