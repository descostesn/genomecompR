.checkGFF <- function(filepath) { # nolint

    if (!isTRUE(all.equal(tools::file_ext(filepath), "gff")))
        stop("The glc peaks should be in GFF format")
}

.filterAnno <- function(currenttable, currentname) { # nolint

    message("\t Filtering ", currentname)
    if (isTRUE(all.equal(currentname, "gencode"))) {

        currenttable$V9 <- "."
        return(currenttable[which(currenttable$V3 == "gene" & 
                                        currenttable$V1 != "chrM"), ])
    }

    idx <- c(grep("Un", currenttable$V1), grep("random", currenttable$V1),
            grep("chrM", currenttable$V1))

    if (!isTRUE(all.equal(length(idx), 0)))
        return(currenttable[-idx, ])
    else
        return(currenttable)
}


#' Create Genomic Compartments Object
#'
#' @usage
#' genomeCompart(peakspathvec, geneannovec)
#'
#' @description
#' This function is the constructor of the genomicCompartments object. It is
#' called by the function buildIntervalsObject.
#'
#' @param peakspathvec Character vector containing the path to the peaks file.
#' @param geneannovec Character vector containing the paths to gene annotation
#' files.
#'
#' @return Returns an object of class `genomicCompartments` containing the
#' processed peaks and gene annotations.
#'
#' @details
#' peakspathvec is made of the GFF file containing glc peaks and of the GFF
#' files containing the ChIP-seq and ATAC-seq peak coordinates. The geneannovec
#' contains path to the gencode, refgene, and refseq gff files.
#'
#' @seealso
#' buildIntervalsObject
#'
#' @examples
#' \dontrun{
#' peakspathquery <- "/path/to/glc.gff"
#' peakspathcategoriesvec <- c(H3K27ac = "/path/to/H3K27ac.gff",
#'        H3K4me1 = "/path/to/H3K4me1.gff", H3K27me3 = "/path/to/H3K27me3.gff",
#'        H3K4me3 = "/path/to/H3K4me3.gff", Suz12 = NA, RING1B = NA,
#'        H3K9me3 = "/path/to/H3K9me3.gff", Ser5P = "/path/to/Ser5P.gff",
#'        Ser2P = "/path/to/Ser2P.gff", ATACSeq = "/path/to/ATACseq.gff")
#' peakspathvec <- c(peakspathquery, peakspathcategoriesvec)
#' geneannovec <- c(gencode = "/path/to/gencode.gff",
#'        refgene = "/path/to/refGeneUCSC.gff",
#'        refseq = "/path/to/refseq.gff")
#' 
#' genomic_compartments <- genomeCompart(peakspathvec, geneannovec)
#' }
#'
genomeCompart <- function(peakspathvec, geneannovec) { # nolint


    ## Reading the glc peaks and converting to GRanges
    message("Reading the glc peaks and converting to GRanges")
    .checkGFF(peakspathvec[1])
    grglc <- buildGR(peakspathvec[1]) # nolint

    ## Reading and filtering annotations:
    ## 1) Gencode: Removing non-gene annotations and chromosome M
    ## 2) refgene: Removing non-canonical chromosomes
    ## 3) refseq: Removing non-canonical chromosome
    message("Reading and filtering genes annotations")
    invisible(sapply(geneannovec, .checkGFF))
    geneannolist <- lapply(geneannovec, read.table)
    geneannolist <- mapply(.filterAnno, geneannolist, names(geneannolist),
            SIMPLIFY = FALSE)

    ## Transforming to a GRangesList
    message("Transforming to a GRangesList")
    geneannolistgr <- lapply(geneannolist, buildGRNoReading) # nolint
    grlist <- GenomicRanges::GRangesList(geneannolistgr)

    new("genomicCompartments",
            referencePeaks = grglc,
            geneAnnotations = grlist)
}
