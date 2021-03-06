% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{makeSequenceGraph}
\alias{makeSequenceGraph}
\title{Given a description of the variants in a mixture, will generate a graph for use with later analysis.}
\usage{
makeSequenceGraph(refRS, diffAllele, position, etype)
}
\arguments{
\item{refRS}{The reference sequence.}

\item{diffAllele}{The variant alleles for each difference.}

\item{position}{The positions in the reference: one based. Multiple events can be specified at a location.}

\item{etype}{The type of each event. Zero for replacement, one for insert, two for delete.}
}
\value{
An opaque object for use with other methods in this package.
}
\description{
Given a description of the variants in a mixture, will generate a graph for use with later analysis.
}
