% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-plot.R
\name{boxplotGlcnacLevels,genomicCompartments-method}
\alias{boxplotGlcnacLevels,genomicCompartments-method}
\title{Generate Boxplot of Glcnac Levels}
\usage{
boxplotGlcnacLevels(theobject, bwpath, outputfolder, includerepeats,
           plotviolin = FALSE)
}
\arguments{
\item{theobject}{An object of class `genomicCompartments`.}

\item{bwpath}{Character string specifying the path to the BigWig file
containing Glcnac levels.}

\item{outputfolder}{Character string specifying the output folder path for
saving results.}

\item{includerepeats}{Logical indicating whether to include repeats in the
analysis.}

\item{plotviolin}{Logical indicating whether to plot to current device
(default is FALSE).}
}
\value{
Save the boxplot to outputfolder.
}
\description{
This method generates boxplots of Glcnac levels across the different genomic
compartments stored in the given genomicCompartments object.
}
\examples{
\dontrun{
# Create a genomicCompartments object
gc_obj <- genomeCompart(peakspathvec, geneannovec)
# Generate boxplots of Glcnac levels
gc_obj <- boxplotGlcnacLevels(gc_obj, bwpath = "path/to/glcnac_levels.bw",
                              outputfolder = "path/to/output_folder",
                              includerepeats = FALSE, plotviolin = FALSE)
}

}
