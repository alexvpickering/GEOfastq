% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srp_utils.R
\name{extract_gsms}
\alias{extract_gsms}
\title{Extract GSMs needed to download RNA-seq data for a series}
\usage{
extract_gsms(gse_text)
}
\arguments{
\item{gse_text}{GSE text returned from \code{\link{crawl_gse}}}
}
\value{
Character vector of sample GSMs for the series \code{gse_name}
}
\description{
Extract GSMs needed to download RNA-seq data for a series
}
\examples{

gse_text <- crawl_gse('GSE111459')
gsm_names <- extract_gsms(gse_text)

}
