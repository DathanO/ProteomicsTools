% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{fisher_dotplot}
\alias{fisher_dotplot}
\title{Title}
\usage{
fisher_dotplot(ftest, path = "examples/")
}
\arguments{
\item{ftest}{The output of fishertest()}

\item{path}{will define directory to store output. by default, files are saved within wd}
}
\value{
A plot of the fisher test enrichment
}
\description{
Title
}
\examples{
data("ftest")
fisher_dotplot(ftest)
}
