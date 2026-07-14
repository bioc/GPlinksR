\name{gp_example_links}
\docType{data}
\alias{gp_example_links}
\title{Example Link Table for the Packaged GPlinksR Demo Subset}
\format{
A data.frame with 7 rows and 3 columns:
\describe{
\item{\code{Peak}}{Peak coordinate in \code{"chr:start-end"} format.}
\item{\code{Gene}}{Gene symbol linked to the peak.}
\item{\code{Src}}{Link source, one of \code{"enh"}, \code{"prom"}, or
\code{"clo"}.}
}
}
\description{
A reference gene-peak link table corresponding to the leading entries in
\code{gp_example_inputs}. The table contains the expected enhancer-, promoter-,
and closest-based links from the package's executable offline example.
}
\usage{data(gp_example_links)}
\examples{
data(gp_example_links)
head(gp_example_links)
table(gp_example_links$Src)
}
