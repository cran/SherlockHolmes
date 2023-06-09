% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sherlock.R
\name{readTitles}
\alias{readTitles}
\title{readTitles}
\usage{
readTitles(titles)
}
\arguments{
\item{titles}{is a character string containing the full path name for a text file containing the titles of the stories in the same order that thney appear in the texts file}
}
\value{
a character vector of titles
}
\description{
read and edit titles to remove blank lines and white space
}
\examples{
titles<-system.file("extdata/contents3.txt",package="SherlockHolmes")
rt<-readTitles(titles)
}
