#' @title uchimeout_help
#'
#' @description Description of 'uchimeout' data.frame returned by the
#' [[rchime()]] function.
#'
#' @export
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
uchimeout_help <- function() {
  message <- paste0(
    "Output description from uchime's manual:", "\n",
    "1	Score	Higher score means more strongly chimeric alignment.", "\n",
    "2	Q	Query Name.", "\n",
    "3	A	Parent A Name.", "\n",
    "4	B	Parent B Name.", "\n",
    "5	T	Top parent (T) Name. This is the closest reference sequence; usually
          either A or B.", "\n",
    "6	IdQM	Percent identity of query and the model (M) constructed as a
              segment of A and a segment of B.", "\n",
    "7	IdQA	Percent identity of Q and A.", "\n",
    "8	IdQA	Percent identity of Q and B.", "\n",
    "9	IdAB	Percent identity of A and B", "\n",
    "10	IdQT	Percent identity of Q and T.", "\n",
    "11	LY	Yes votes in left segment.", "\n",
    "12	LN	No votes in left segment.", "\n",
    "13	LA	Abstain votes in left segment.", "\n",
    "14	RY	Yes votes in right segment.", "\n",
    "15	RN	No votes in right segment.", "\n",
    "16	RA	Abstain votes in right segment.", "\n",
    "17	Div	Divergence, defined as (IdQM - IdQT).", "\n",
    "18	YN	Y, N or ?, indicating whether the query was classified as", "\n",
    "chimeric (Y), not chimeric (N) or borderline case (?).", "\n",
    "The query is classified as chimeric if h >= threshold specified by", "\n",
    "minh, Div > minimum divergence specified by mindiv and the number", "\n",
    "of diffs ( (Y+N+A) in each segment (L and R) is greater than the", "\n",
    " minimum specified by -mindiffs. A query is unclassified if the", "\n",
    "maxh > h > minh, i.e. maxh is the maximum score for a non-chimera,", "\n",
    "and minh is the minimum score for a chimera; in between is", "\n",
    "unclassified.", "\n",
    collapse = ""
  )
  cat(message)
}
