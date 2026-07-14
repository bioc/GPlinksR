#' Example Link Table for the Packaged GPlinksR Demo Subset
#'
#' A reference gene-peak link table corresponding to the leading entries in
#' `gp_example_inputs`. The table contains the expected enhancer-, promoter-,
#' and closest-based links from the package's executable offline example.
#'
#' @format A data.frame with 7 rows and 3 columns:
#' \describe{
#'   \item{Peak}{Peak coordinate in `"chr:start-end"` format.}
#'   \item{Gene}{Gene symbol linked to the peak.}
#'   \item{Src}{Link source, one of `"enh"`, `"prom"`, or `"clo"`.}
#' }
"gp_example_links"
