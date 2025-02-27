% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{buildGRNoReading}
\alias{buildGRNoReading}
\title{Build GRanges Object from Data Frame}
\usage{
buildGRNoReading(fi)
}
\arguments{
\item{fi}{A data frame with columns representing genomic coordinates and
other relevant information in GFF format.
The expected columns are:
\itemize{
  \item \code{V1} - Chromosome names
  \item \code{V4} - Start positions
  \item \code{V5} - End positions
  \item \code{V9} - Feature names
  \item \code{V7} - Strand information
}}
}
\value{
A `GRanges` object with the specified genomic coordinates and
feature names.
}
\description{
This function constructs a `GRanges` object from a data frame.
}
\examples{
\dontrun{
# Example data frame
df <- data.frame(V1 = "chr1", V4 = 1000, V5 = 2000, V9 = "gene1", V7 = "+")
# Build GRanges object
gr <- buildGRNoReading(df)
}

}
