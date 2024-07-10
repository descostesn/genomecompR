#' Build GRanges Object from Data Frame
#'
#' @description
#' This function constructs a `GRanges` object from a data frame.
#'
#' @usage
#' buildGRNoReading(fi)
#'
#' @param fi A data frame with columns representing genomic coordinates and
#' other relevant information.
#' The expected columns are:
#' \itemize{
#'   \item \code{V1} - Chromosome names
#'   \item \code{V4} - Start positions
#'   \item \code{V5} - End positions
#'   \item \code{V9} - Feature names
#'   \item \code{V7} - Strand information
#' }
#'
#' @return A `GRanges` object with the specified genomic coordinates and
#' feature names.
#'
#' @examples
#' \dontrun{
#' # Example data frame
#' df <- data.frame(V1 = "chr1", V4 = 1000, V5 = 2000, V9 = "gene1", V7 = "+")
#' # Build GRanges object
#' gr <- buildGRNoReading(df)
#' }
#'
buildGRNoReading <- function(fi) { # nolint

    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = fi$V9), strand = fi$V7)
    return(gr)
}

#' Build GRanges Object from File
#'
#' @description
#' This function reads a file containing genomic coordinates and constructs a
#' `GRanges` object.
#'
#' @usage
#' buildGR(currentpath)
#'
#' @param currentpath A character string representing the GFF file path to read.
#' The file should contain tabular data with the following columns:
#' \itemize{
#'   \item \code{V1} - Chromosome names
#'   \item \code{V4} - Start positions
#'   \item \code{V5} - End positions
#'   \item \code{V9} - Feature names
#'   \item \code{V7} - Strand information
#' }
#'
#' @return A `GRanges` object with the specified genomic coordinates and
#' feature names.
#'
#' @examples
#' \dontrun{
#' # Example file path
#' path <- "path/to/genomic/coordinates/file.gff"
#' # Build GRanges object
#' gr <- buildGR(path)
#' }
#'
buildGR <- function(currentpath) { # nolint

    message("\t Processing ", currentpath)
    fi <- read.table(currentpath, stringsAsFactors = FALSE)
    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = fi$V9), strand = fi$V7)
    return(gr)
}

.filterChrom <- function(currenttable) { # nolint

    idx <- c(grep("Un", currenttable$V1), grep("random", currenttable$V1),
            grep("chrM", currenttable$V1))

    if (!isTRUE(all.equal(length(idx), 0)))
        return(currenttable[-idx, ])
    else
        return(currenttable)
}

.overlapGlucnacComp <- function(theobject, includerepeats) { # nolint

    ## Make a GRangesList with all the compartments
    message("Preparing data.")
    complist <- aslist(theobject, includerepeats) # nolint
    compgrlist <- GenomicRanges::GRangesList(complist)

    ## Retrieve the GR of glucnac peaks
    querygr <- getRefPeaks(theobject) # nolint

    ## Perform the overlap
    message("Performing overlap")
    resultoverlap <- GenomicRanges::findOverlaps(querygr, compgrlist,
        ignore.strand = TRUE)

    return(resultoverlap)
}

#' Validate Parameters for Genomic Analysis
#'
#' @description
#' This function validates various parameters required for genomic analysis to
#' ensure they meet the expected criteria.
#'
#' @usage
#' checkParams <- function(peakspathqueryvec, glcnacbwvec, querynamevec,
#'        geneannovec, peakspathcategoriesvec, repeatsannovec, countstable,
#'        includerepeats, countsannotype)
#'
#' @param peakspathqueryvec A vector of file paths to peak files for the
#' queries.
#' @param glcnacbwvec A vector of file paths to glucnac bigwig files.
#' @param querynamevec A vector of names corresponding to the query files.
#' @param geneannovec A vector of file paths to gene annotation files,
#' expected to be in the order of gencode, refGene, and refseq.
#' @param peakspathcategoriesvec A vector of names for each peak file, expected
#' to contain strings identifying categories such as H3K27ac, H3K4me1, etc.
#' @param repeatsannovec A vector of file paths to repeat annotation files,
#' expected to contain strings identifying LINE, LTR, and SINE.
#' @param countstable A vector containing the path to a counts table file.
#' @param includerepeats A logical value indicating whether to include repeat
#' annotations.
#' @param countsannotype A string specifying the annotation type for the counts
#' table, either "ensembl" or "entrez".
#'
#' @return None. The function stops with an error message if any validation
#' check fails.
#'
#' @details This function performs the following validation checks:
#' \itemize{
#'   \item Ensures there are no more than two experiments.
#'   \item Checks that the number of peak files matches the number of bigwig
#' files.
#'   \item Verifies that the number of query names matches the number of files.
#'   \item Validates that gene annotation files are gencode, refGene, and
#' refseq.
#'   \item Confirms that category peak files contain the expected strings for
#' categories like H3K27ac, H3K4me1, etc.
#'   \item If repeat annotations are included, checks that the repeat files
#' contain the expected strings for LINE, LTR, and SINE.
#'   \item Ensures only one counts table is provided.
#'   \item Validates that the counts table annotation type is either "ensembl"
#' or "entrez".
#' }
#'
#' @examples
#' \dontrun{
#' peakspathqueryvec <- c("path/to/peak1", "path/to/peak2")
#' glcnacbwvec <- c("path/to/bw1", "path/to/bw2")
#' querynamevec <- c("query1", "query2")
#' geneannovec <- c("gencode", "refGene", "refseq")
#' peakspathcategoriesvec <- c("H3K27ac_path", "H3K4me1_path", "H3K27me3_path",
#' "H3K4me3_path", "Suz12_path", "RING1B_path", "H3K9me3_path", "Ser5P_path",
#' "Ser2P_path", "ATAC_path")
#' repeatsannovec <- c("LINE_path", "LTR_path", "SINE_path")
#' countstable <- c("path/to/counts_table")
#' includerepeats <- TRUE
#' countsannotype <- "ensembl"
#' checkParams(peakspathqueryvec, glcnacbwvec, querynamevec, geneannovec,
#' peakspathcategoriesvec, repeatsannovec, countstable, includerepeats,
#' countsannotype)
#' }
#'
checkParams <- function(peakspathqueryvec, glcnacbwvec, querynamevec, # nolint
        geneannovec, peakspathcategoriesvec, repeatsannovec, countstable,
        includerepeats, countsannotype) {

    if (length(peakspathqueryvec) > 2)
        stop("The script is currently designed to handle two experiments.")

    if (!isTRUE(all.equal(length(peakspathqueryvec), length(glcnacbwvec))))
        stop("The number of peak files and bigwig files does not match.")

    if (!isTRUE(all.equal(length(peakspathqueryvec), length(querynamevec))))
        stop("The number of names of query does not match the number of files.")

    if (!(grepl("gencode", geneannovec[1]) && grepl("refGene", geneannovec[2])
                && grepl("refseq", geneannovec[3])))
        stop("Verify that your gene annotations are gencode, refGene, and ",
                "refseq. Each file should contain these strings.")

    if (grepl("H3K27ac", peakspathcategoriesvec[2]) &&
            grepl("H3K4me1", peakspathcategoriesvec[3]) &&
            grepl("H3K27me3", peakspathcategoriesvec[4]) &&
            grepl("H3K4me3", peakspathcategoriesvec[5]) &&
            grepl("Suz12", peakspathcategoriesvec[6]) &&
            grepl("RING1B", peakspathcategoriesvec[7]) &&
            grepl("H3K9me3", peakspathcategoriesvec[8]) &&
            grepl("Ser5P", peakspathcategoriesvec[9]) &&
            grepl("Ser2P", peakspathcategoriesvec[10]) &&
            grepl("ATAC", peakspathcategoriesvec[11]))
        stop("Verify that your category peaks are H3K27ac, H3K4me1, H3K27me3,",
                " H3K4me3, Suz12, RING1B, H3K9me3, Ser5P, Ser2P, ATACSeq. ",
                "Each file should contain these strings.")

    if (includerepeats)
        if (!(grepl("LINE", repeatsannovec[1]) &&
                    grepl("LTR", repeatsannovec[2]) &&
                    grepl("SINE", repeatsannovec[3])))
            stop("Verify that your repeats are LINE, LTR, and SINE. Each ",
                    "file should contain these strings.")

    if (!isTRUE(all.equal(length(countstable), 1)))
        stop("Only one counts table should be provided.")

    if (!isTRUE(all.equal(countsannotype, "ensembl")) &&
        !isTRUE(all.equal(countsannotype, "entrez")))
        stop("The first column of the count table file should be ensembl ",
            " or entrez IDs")
}

#' Filter Chromosomes and Strand Information in GRanges
#'
#' @description
#' This function filters out unwanted chromosomes and standardizes strand
#' information in a `GRanges` object.
#'
#' @usage
#' filterChromAndStrand(currentgr)
#'
#' @param currentgr A `GRanges` object containing genomic ranges to be filtered.
#'
#' @return A filtered `GRanges` object with unwanted chromosomes removed and
#' strand information standardized.
#'
#' @details This function performs the following operations:
#' \itemize{
#'   \item If the strand information contains neither "+" and "-" strands, all
#' strand information is set to "+".
#'   \item Chromosomes with names containing "Un", "random", or "chrM" are
#' removed.
#'   \item The sequence levels of the resulting `GRanges` object are updated to
#' include only the remaining chromosomes.
#' }
#'
#' @examples
#' \dontrun{
#' # Example GRanges object
#' gr <- GRanges(seqnames = c("chr1", "chr2", "chrUn", "chrM"),
#'               ranges = IRanges(start = c(1, 2, 3, 4), end = c(5, 6, 7, 8)),
#'               strand = c("+", "-", "+", "*"))
#' # Filter the GRanges object
#' filtered_gr <- filterChromAndStrand(gr)
#' }
#'
filterChromAndStrand <- function(currentgr) { # nolint

    if (!isTRUE(all.equal(unique(
        as.character(BiocGenerics::strand(currentgr))), "+")) &&
            !isTRUE(all.equal(unique(
                as.character(BiocGenerics::strand(currentgr))), "-")))
        BiocGenerics::strand(currentgr) <- "+"

    chromvec <- as.character(GenomeInfoDb::seqnames(currentgr))
    idx <- c(grep("Un", chromvec), grep("random", chromvec),
            grep("chrM", chromvec))

    if (!isTRUE(all.equal(length(idx), 0))) {

        currentgr <- currentgr[-idx, ]
        GenomeInfoDb::seqlevels(currentgr) <- unique(
            as.character(GenomeInfoDb::seqnames(currentgr)))
        return(currentgr)
    } else {
        return(currentgr)
    }
}

#' Build Genomic Intervals Object
#'
#' @description
#' This function constructs a `genomeCompart` object, defines various genomic
#' features, and optionally includes repeat annotations.
#'
#' @usage
#' buildIntervalsObject(peakspathvec, geneannos, includerepeats)
#'
#' @param peakspathvec A named vector of file paths to peak files for various
#' ChIP-seq and ATAC-seq data. Expected names include: "H3K27ac", "Ser5P",
#' "Ser2P", "H3K4me3", "H3K27me3", "ATACSeq", "Suz12", "RING1B", "H3K9me3".
#' @param geneannos A list of gene annotation files to be used for defining
#' non-redundant gene sets and other features.
#' @param includerepeats A logical value indicating whether to include repeat
#' annotations (LINE, LTR, SINE).
#'
#' @return An updated `genomeCompart` object with defined genomic features and,
#' optionally, repeat annotations.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Builds the initial `genomeCompart` object.
#'   \item Defines a non-redundant set of gencode genes by removing overlapping
#' genes, genes closer than 1kb to each other, and genes shorter than 2kb.
#'   \item Identifies active promoters using H3K27ac at TSS±1Kb.
#'   \item Identifies promoters with Ser5P (initiation) at TSS±1Kb.
#'   \item Identifies elongating genes by overlapping Ser2P from TSS+1kb to TES.
#'   \item Defines termination sites at TES+1kb.
#'   \item Identifies bivalent promoters using H3K4me3 and H3K27me3 at TSS±1kb.
#'   \item Defines active enhancers using H3K27ac, H3K4me1, and ATAC-seq,
#' ensuring no overlap with gene annotations.
#'   \item Defines poised enhancers using H3K27me3, H3K4me1, and Suz12,
#' ensuring no overlap with gene annotations.
#'   \item Defines polycomb domains using Suz12 and RING1B overlaps.
#'   \item Defines heterochromatin using H3K9me3.
#'   \item Optionally includes repeat annotations (SINE, LINE, LTR) if
#' `includerepeats` is TRUE.
#' }
#'
#' @examples
#' \dontrun{
#' peakspathvec <- c(H3K27ac = "path/to/H3K27ac.gff",
#'      Ser5P = "path/to/Ser5P.gff", Ser2P = "path/to/Ser2P.gff",
#'      H3K4me3 = "path/to/H3K4me3.gff", H3K27me3 = "path/to/H3K27me3.gff",
#'      ATACSeq = "path/to/ATACSeq.gff", Suz12 = "path/to/Suz12.gff",
#'      RING1B = "path/to/RING1B.gff", H3K9me3 = "path/to/H3K9me3.gff")
#' geneannos <- list("path/to/gencode.gff", "path/to/refGene.gff",
#' "path/to/refseq.gff")
#' includerepeats <- FALSE
#' gc <- buildIntervalsObject(peakspathvec, geneannos, includerepeats)
#' }
#'
buildIntervalsObject <- function(peakspathvec, geneannos, includerepeats, # nolint
    repeatsannovec = NA) {

    ## Verify that repeatsannovec is provided if includerepeats is TRUE
    if (includerepeats && is.na(repeatsannovec))
        stop("repeatsannovec should be provided if includerepeats is TRUE.")

    ## Build the initial object
    gc <- genomeCompart(peakspathvec, geneannos) # nolint

    ## Define a non redundant set of gencode genes: No overlap, remove genes
    ## closer to 1kb from one another, remove genes that are too short i.e.
    ## length < 2kb.
    gc <- annoCuration(gc, geneannos) # nolint

    ## Using the non redundant set of genes, define active promoter as K27ac in
    ## TSS-/+1Kb
    if (!is.na(peakspathvec["H3K27ac"]))
      gc <- overlapOnGenes(gc, "TSS", peaklabel = "H3K27ac", # nolint
            peakpath = peakspathvec["H3K27ac"])
    else
      message("Active promoters could not be defined, H3K27ac bigwig not found")

    ## Using the non redundant set of genes, define promoter with Ser5P
    ## (initiation) in TSS-/+1Kb
    if (!is.na(peakspathvec["Ser5P"]))
      gc <- overlapOnGenes(gc, "TSS", peaklabel = "Ser5P",  # nolint
            peakpath = peakspathvec["Ser5P"])
    else
      message("Promoters could not be defined, Ser5P bigwig not found")

    ## Using the non redundant set of genes, define elongating genes as TSS+1kb
    ## to TES overlapping Ser2P.
    if (!is.na(peakspathvec["Ser2P"]))
      gc <- overlapOnGenes(gc, "GBTES", peaklabel = "Ser2P", # nolint
            peakpath = peakspathvec["Ser2P"])
    else
      message("Elongation could not be defined, Ser2P bigwig not found")

    ## Using the non redundant set of genes, define termination sites as
    ## TES+1kb
    gc <- overlapOnGenes(gc, "TES") # nolint

    ## Using the non redundant set of genes, define bivalent promoters: H3K4me3
    ## and H3K27me3 at TSS+/- 1kb
    if (!is.na(peakspathvec["H3K4me3"]) && !is.na(peakspathvec["H3K27me3"]))
      gc <- bivalentPromoters(gc, peakspathvec) # nolint
    else
      message("Bivalent promoters could not be defined, H3K27me3 and H3K4me3",
        " bigwigs not found")

    ## Define active enhancer: H3K27ac/H3K4me1/ATAC-seq. They should not
    ## overlap with the combination of UCSC refGene, NCBI RefSeq, and GENCODE
    ## VM25 that are in the geneAnnotations slot.
    if (!is.na(peakspathvec["H3K27ac"]) && !is.na(peakspathvec["H3K4me1"]) &&
      !is.na(peakspathvec["ATACSeq"]))
      gc <- enhancer(gc, peakspathvec["H3K27ac"], peakspathvec["H3K4me1"], # nolint
            peakspathvec["ATACSeq"], "H3K27ac", "H3K4me1", "ATACSeq", "active")
    else
      message("Active enhancers could not be defined, need H3K27ac, H3K4me1",
        " and ATAC-seq bigwigs")

    ## Define poised enhancer: H3K27me3/H3K4me1/PRC2. They should not overlap
    ## with the combination of UCSC refGene, NCBI RefSeq, and GENCODE VM25.
    if (!is.na(peakspathvec["H3K27me3"]) && !is.na(peakspathvec["H3K4me1"]) &&
      !is.na(peakspathvec["Suz12"]))
      gc <- enhancer(gc, peakspathvec["H3K27me3"], peakspathvec["H3K4me1"], # nolint
            peakspathvec["Suz12"], "H3K27me3", "H3K4me1", "Suz12", "poised")
    else
      message("Poised enhancers could not be defined, need H3K27me3, H3K4me1",
        " and Suz12 bigwigs")


    ## Define polycomb domain: Suz12 and RING1B overlap
    if (!is.na(peakspathvec["Suz12"]) && !is.na(peakspathvec["RING1B"]))
      gc <- polycombsDomains(gc, peakspathvec["Suz12"], peakspathvec["RING1B"], # nolint
            "Suz12", "RING1B", "polycombs")
    else
      message("Polycomb domains could not be defined, need Suz12 and Ring1b",
        " bigwigs")

    ## Define heterochromatin: H3K9me3
    if (!is.na(peakspathvec["H3K9me3"]))
      setHeterochromatin(gc) <- filterChromAndStrand(buildGR( # nolint
                    peakspathvec["H3K9me3"]))
    else
      message("Heterochromatin coordinates could not be defined, need H3K9me3",
        " bigwig")

    ## Define SINE, LINE and LTR
    if (includerepeats) {
        setSINE(gc) <- filterChromAndStrand(buildGR(repeatsannovec["SINE"])) # nolint
        setLINE(gc) <- filterChromAndStrand(buildGR(repeatsannovec["LINE"])) # nolint
        setLTR(gc) <- filterChromAndStrand(buildGR(repeatsannovec["LTR"])) # nolint
    }
    return(gc)
}

createfolder <- function(outfold) {
    if (!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
}

#' Retrieve Gene Conversion Table
#'
#' @description
#' This function retrieves a conversion table of gene identifiers from
#' Ensembl BioMart.
#'
#' @param biomartstr A string specifying the BioMart database to use (e.g.,
#' "ENSEMBL_MART_ENSEMBL").
#' @param datasetstr A string specifying the dataset to use within the BioMart
#' database (e.g., "hsapiens_gene_ensembl").
#' @param hoststr A string specifying the host URL for BioMart (e.g.,
#' "http://www.ensembl.org").
#' @param alternativemirroropt A string specifying an alternative mirror option
#' for BioMart, if needed.
#'
#' @return A data frame containing the conversion table with columns for RefSeq,
#' Entrez, and Ensembl gene identifiers.
#'
#' @details The function filters the retrieved data to include only rows with
#' non-empty and non-NA gene IDs.
#'
#' @examples
#' \dontrun{
#' biomartstr <- "ENSEMBL_MART_ENSEMBL"
#' datasetstr <- "hsapiens_gene_ensembl"
#' hoststr <- "http://www.ensembl.org"
#' alternativemirroropt <- "useast"
#' geneconvert <- retrieveconversiontab(biomartstr, datasetstr, hoststr,
#' alternativemirroropt)
#' head(geneconvert)
#' }
#'
retrieveconversiontab <- function(biomartstr, datasetstr, hoststr,
  alternativemirroropt) {

    ensembl <- tryUseMart(biomart = biomartstr, # nolint
            dataset = datasetstr,
            host = hoststr,
            alternativemirror = alternativemirroropt)
    message("\t Retrieving gene info")
    infovec <- c("refseq_mrna", "refseq_ncrna", "entrezgene_id",
        "ensembl_gene_id")
    geneconvert <- tryGetBM(infovec, ensembl) # nolint
    colnames(geneconvert) <- c("", "", "", "")
    idxmrna <- which(geneconvert[, 1] != "" & !is.na(geneconvert[, 3]))
    idxnc <- which(geneconvert[, 2] != "" & !is.na(geneconvert[, 3]))
    geneconvert <- as.data.frame(rbind(geneconvert[idxmrna, c(1, 3, 4)],
                    geneconvert[idxnc, c(2, 3, 4)]))
    colnames(geneconvert) <- c("refseq", "entrez", "ensembl")
    return(geneconvert)
}
