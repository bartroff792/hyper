% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperMCI.R
\name{AMO}
\alias{AMO}
\title{AMO intervals}
\usage{
AMO(N, n, a)
}
\arguments{
\item{N}{population size}

\item{n}{sample size}

\item{a}{alpha}
}
\value{
AMO intervals
}
\description{
Create alpha max optimal acceptance intervals
}
\examples{
AMO_interval <- AMO(100,10,0.05);
}
