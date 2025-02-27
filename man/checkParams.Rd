% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{checkParams}
\alias{checkParams}
\title{Validate Parameters given as input to different functions of the package}
\usage{
checkParams(peakspathqueryvec, glcnacbwvec, querynamevec,
       geneannovec, peakspathcategoriesvec, repeatsannovec, countstable,
       includerepeats, countsannotype)
}
\arguments{
\item{peakspathqueryvec}{A vector of file paths to peak files in GFF format.}

\item{glcnacbwvec}{A vector of file paths to glc bigwig files.}

\item{querynamevec}{A vector of names corresponding to the peak files.}

\item{geneannovec}{A vector of file paths to gene annotation files,
expected to be in the order of gencode, refGene, and refseq.}

\item{peakspathcategoriesvec}{A vector of names for each peak file, expected
to contain strings identifying categories such as H3K27ac, H3K4me1, etc.}

\item{repeatsannovec}{A vector of file paths to repeat annotation files in
GFF format. The file names are expected to contain the substrings LINE,
LTR, and SINE.}

\item{countstable}{A vector containing the path to a count table file.}

\item{includerepeats}{A logical value indicating whether to include repeat
annotations.}

\item{countsannotype}{A string specifying the annotation type for the counts
table, either "ensembl" or "entrez".}
}
\value{
None. The function stops with an error message if any validation
check fails.
}
\description{
This function validates various parameters required by several functions to
ensure they meet the expected criteria.
}
\details{
This function performs the following validation checks:
\itemize{
  \item Ensures there are no more than two files in peakspathqueryvec.
  \item Checks that the number of peak files matches the number of bigwig
files.
  \item Verifies that the number of query names matches the number of peak
files.
  \item Validates that gene annotation files are ordered as gencode, refGene,
and refseq.
  \item Confirms that category peak files contain the expected strings for
categories like H3K27ac, H3K4me1, etc.
  \item If repeat annotations are included, checks that the repeat files
contain the expected substrings LINE, LTR, and SINE.
  \item Ensures only one counts table is provided.
  \item Validates that the counts table annotation type is either "ensembl"
or "entrez".
}
}
\examples{
\dontrun{
peakspathqueryvec <- c("path/to/peak1.gff", "path/to/peak2.gff")
glcnacbwvec <- c("path/to/bw1", "path/to/bw2")
querynamevec <- c("query1", "query2")
geneannovec <- c("gencode", "refGene", "refseq")
peakspathcategoriesvec <- c("H3K27ac_path", "H3K4me1_path", "H3K27me3_path",
"H3K4me3_path", "Suz12_path", "RING1B_path", "H3K9me3_path", "Ser5P_path",
"Ser2P_path", "ATAC_path")
repeatsannovec <- c("LINE_path", "LTR_path", "SINE_path")
countstable <- c("path/to/counts_table")
includerepeats <- TRUE
countsannotype <- "ensembl"
checkParams(peakspathqueryvec, glcnacbwvec, querynamevec, geneannovec,
peakspathcategoriesvec, repeatsannovec, countstable, includerepeats,
countsannotype)
}

}
