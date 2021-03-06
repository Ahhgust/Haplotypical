% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fastCloseIndelPair}
\alias{fastCloseIndelPair}
\title{Returns the size of indel between 2 reads. Only single indels are supported}
\usage{
fastCloseIndelPair(query, target, ignoreHomopolymers = FALSE)
}
\arguments{
\item{query}{(query sequence)}

\item{target}{(target sequence)}

\item{ignoreHomopolymers}{(consecutive sequences of the same letter is treated as just 1 instance of that letter)}
}
\description{
This function computes a naive (and fast)
"indel" distance. i.e., if the two reads are separated by a single indel whose
length is the distance between the two strings.
This function is an optimization of the common case of sequence comparisons
between reads from massively parallel sequencing that are from the same PCR/sequencing assay
in which case the distance between strings is almost always 1.
It optionally ignores homopolymers (e.g., AAAT == AT == distance of 0 between)
The size of the indel is returned iff they are separated by an indel.
otherwise -1 is returned
}
