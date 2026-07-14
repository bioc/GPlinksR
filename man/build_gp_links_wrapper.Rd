\name{build_gp_links_wrapper}
\alias{build_gp_links_wrapper}
\title{Build Gene-Peak Links from MAE or SCE Input}
\usage{
build_gp_links_wrapper(
    x,
    peak_experiment = NULL,
    gene_experiment = NULL,
    peak_col = NULL,
    gene_col = NULL,
    ver = 19,
    enh_file = NULL,
    enh_coords_file = NULL,
    gene_map = NULL,
    mart = NULL
)
}
\arguments{
\item{x}{A \code{MultiAssayExperiment}, \code{SingleCellExperiment},
\code{SummarizedExperiment}, or \code{RangedSummarizedExperiment}.}

\item{peak_experiment}{For \code{MultiAssayExperiment}, the experiment name
holding peak features. For \code{SingleCellExperiment}, the \code{altExp()}
name holding peak features.}

\item{gene_experiment}{For \code{MultiAssayExperiment}, the experiment name
holding gene features. Ignored for \code{SingleCellExperiment}, where the main
experiment is used for genes.}

\item{peak_col}{Optional \code{rowData()} column containing peak coordinates in
\code{"chr:start-end"} format.}

\item{gene_col}{Optional \code{rowData()} column containing gene symbols.}

\item{ver}{Integer (17, 18, or 19). PEREGRINE enhancer-gene link version.}

\item{enh_file}{Optional local path to enhancer-gene link file (\code{.tsv}).}

\item{enh_coords_file}{Optional local path to the corresponding PEREGRINE
enhancer-coordinate file passed to \code{build_gp_links()}.}

\item{gene_map}{Optional data.frame with columns \code{hgnc_symbol} and
\code{hgnc_id} passed to \code{build_gp_links()}.}

\item{mart}{Optional \code{biomaRt::Mart} object passed to
\code{build_gp_links()}. Supply an object created with
\code{biomaRt::useEnsembl()} to use an Ensembl mirror or custom host.}
}
\value{
A data.frame with columns returned by \code{build_gp_links()}.
}
\description{
Wrapper around \code{build_gp_links()} that extracts peak coordinates and gene
symbols from a \code{MultiAssayExperiment} or \code{SingleCellExperiment},
then runs the standard GPlinksR workflow.
}
\details{
Peak features are extracted from \code{rowRanges()}, \code{rowData()}, or row
names. Gene symbols are extracted from \code{rowData()} or row names.
}
\examples{
if (requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(seq_len(12), nrow = 3)),
        rowData = S4Vectors::DataFrame(symbol = c("TP53", "GATA1", "SPI1"))
    )

    SingleCellExperiment::altExp(sce, "ATAC") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = matrix(seq_len(12), nrow = 3)),
            rowData = S4Vectors::DataFrame(
                PeakRegion = c(
                    "chr1:1000-1200",
                    "chr1:5000-5200",
                    "chr2:9000-9300"
                )
            )
        )

    sce
    # Example wrapper call:
    # build_gp_links_wrapper(
    #     x = sce, peak_experiment = "ATAC", gene_col = "symbol"
    # )
}
}
