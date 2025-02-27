% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{filterChromAndStrand}
\alias{filterChromAndStrand}
\title{Filter Chromosomes and Strand Information in GRanges}
\usage{
filterChromAndStrand(currentgr)
}
\arguments{
\item{currentgr}{A `GRanges` object containing genomic ranges to be filtered.}
}
\value{
A filtered `GRanges` object with unwanted chromosomes removed and
strand information standardized.
}
\description{
This function filters out unwanted chromosomes and standardizes strand
information in a `GRanges` object.
}
\details{
This function performs the following operations:
\itemize{
  \item If the strand information contains neither "+" and "-" strands, all
strand information is set to "+".
  \item Chromosomes with names containing "Un", "random", or "chrM" are
removed.
  \item The sequence levels of the resulting `GRanges` object are updated to
include only the remaining chromosomes.
}
}
\examples{
\dontrun{
# Example GRanges object
gr <- GRanges(seqnames = c("chr1", "chr2", "chrUn", "chrM"),
              ranges = IRanges(start = c(1, 2, 3, 4), end = c(5, 6, 7, 8)),
              strand = c("+", "-", "+", "*"))
# Filter the GRanges object
filtered_gr <- filterChromAndStrand(gr)
}

}
