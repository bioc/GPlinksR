#' GPlinksR: Gene-Peak Regulatory Network Construction
#'
#' @description
#' GPlinksR constructs gene-peak regulatory networks for ATAC-RNA integration.
#' It combines enhancer-gene links, promoter overlaps, and closest-gene
#' mappings, and accepts either direct feature vectors or common Bioconductor
#' container objects.
#'
#' @section Main functions:
#' * [build_gp_links()] constructs links from peak coordinates and gene symbols.
#' * [build_gp_links_wrapper()] extracts inputs from Bioconductor containers.
#' * [get_peregrine_file()] downloads and caches PEREGRINE enhancer-gene data.
#'
#' @docType package
#' @name GPlinksR-package
#' @aliases GPlinksR
#' @keywords internal
"_PACKAGE"
