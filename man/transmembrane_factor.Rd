% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/API_functions.R
\name{transmembrane_factor}
\alias{transmembrane_factor}
\title{create a new factor column based on the number of transmembrane domains.}
\usage{
transmembrane_factor(df, colname, newcolname, treshold = 3)
}
\arguments{
\item{df}{with a column containing transmembrane domains}

\item{colname}{string name of the new column created by the function}

\item{newcolname}{string of the new colname}

\item{treshold}{number of transmembrane to seperate. It will create 3 groups, 0, 1-treshold, and >treshold.}
}
\value{
a new dataframe with a new column
}
\description{
create a new factor column based on the number of transmembrane domains.
}
\examples{
transmembrane_factor(datamined, "TMregions", "TMgroups")
}
