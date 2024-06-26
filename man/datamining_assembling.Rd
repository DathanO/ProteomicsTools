% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/API_functions.R
\name{datamining_assembling}
\alias{datamining_assembling}
\title{Incorporate a new column from the datamining functions above}
\usage{
datamining_assembling(df, colname, fun)
}
\arguments{
\item{df}{containing "name" as colname}

\item{colname}{in strings of the datamined column}

\item{fun}{name of the function to use}
}
\value{
df containing the new column
}
\description{
Incorporate a new column from the datamining functions above
}
\examples{
datamining_assembling(head(datamined), "PI", get_nextprot_PI)
}
