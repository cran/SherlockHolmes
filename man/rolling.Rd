% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sherlock.R
\name{rolling}
\alias{rolling}
\title{rolling}
\usage{
rolling(freqs, titles.vec, windowPct = 0.1, odir, verbose)
}
\arguments{
\item{freqs}{return value of frequency()}

\item{titles.vec}{character vector containing the titles of the stories}

\item{windowPct}{a numeric control size of plot window}

\item{odir}{character string containing the full path name for the output directory}

\item{verbose}{Boolean if TRUE print informative or diagnostic messages to console}
}
\value{
returns noo value, but has side effect of generating graphs
}
\description{
compute rolling average of ratio of number of occurrences of query string divided by total number of words
}
\examples{
rol<-rolling(freqs,titles.vec,windowPct=0.10,odir=tempdir(),verbose=FALSE)
}
