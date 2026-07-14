\name{gp_example_inputs}
\docType{data}
\alias{gp_example_inputs}
\title{Example Gene and Peak Inputs for GPlinksR}
\format{
A list with two elements:
\describe{
\item{\code{pk}}{A character vector of 300 peak coordinates in
    \code{"chr:start-end"} format.}
\item{\code{gn}}{A character vector of 100 gene symbols.}
}
}
\description{
A small example dataset for demonstrating \code{GPlinksR} workflows.
The object contains 300 peak coordinates and 100 gene symbols selected from a
KPMP-derived kidney single-cell analysis. The first few entries yield
enhancer-, promoter-, and closest-based links in the executable examples.
}
\usage{data(gp_example_inputs)}
\examples{
data(gp_example_inputs)
length(gp_example_inputs$pk)
length(gp_example_inputs$gn)
}
