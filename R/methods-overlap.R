.defineIntervals <- function(grcuratedanno, regionlabel) { # nolint

    if (isTRUE(all.equal(regionlabel, "TSS")))
        gr <- GenomicRanges::GRanges(
                seqnames = as.character(GenomicRanges::seqnames(grcuratedanno)),
                ranges = IRanges::IRanges(start =
                    GenomicRanges::start(grcuratedanno) - 1000,
                        end = GenomicRanges::start(grcuratedanno) + 1000),
                strand = GenomicRanges::strand(grcuratedanno))
    else if (isTRUE(all.equal(regionlabel, "GBTES")))
        gr <- GenomicRanges::GRanges(
                seqnames = as.character(GenomicRanges::seqnames(grcuratedanno)),
                ranges = IRanges::IRanges(start =
                    GenomicRanges::start(grcuratedanno) + 1000,
                        end = GenomicRanges::end(grcuratedanno)),
                strand = GenomicRanges::strand(grcuratedanno))
    else if (isTRUE(all.equal(regionlabel, "TES")))
        gr <- GenomicRanges::GRanges(
                seqnames = as.character(GenomicRanges::seqnames(grcuratedanno)),
                ranges = IRanges::IRanges(
                    start = GenomicRanges::end(grcuratedanno), 
                    end = GenomicRanges::end(grcuratedanno) + 50),
                strand = GenomicRanges::strand(grcuratedanno))
    else
        stop("The label ", regionlabel, " is not recognized")

    return(gr)
}

#' Overlap Peaks on Genes
#'
#' @description
#' This method performs overlap analysis between peaks and curated annotations
#' at specified genomic regions (TSS or TES). It assigns overlapping peaks to
#' specific compartments based on the provided labels.
#'
#' @usage
#' overlapOnGenes(theobject, regionlabel, peaklabel = NA, peakpath = NA)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param regionlabel Character specifying the genomic region label ("TSS",
#' "GBTES", or "TES").
#' @param peaklabel Character specifying the peak label ("H3K27ac", "Ser5P",
#' or "Ser2P").
#' @param peakpath Path to the file containing peaks data. NA if regionlabel is
#' "TES".
#'
#' @return Returns the modified `genomicCompartments` object with assigned
#' compartments.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Perform overlap analysis on TSS with H3K27ac peaks
#' overlapOnGenes(gc_obj, regionlabel = "TSS", peaklabel = "H3K27ac",
#' peakpath = "path/to/peaks.gff")
#' }
#'
setMethod( # nolint

        f = "overlapOnGenes",

        signature = "genomicCompartments",

        definition = function(theobject, regionlabel, peaklabel = NA,
            peakpath = NA) {

            ## Check the Object
            validObject(theobject)

            ## Retrieve the curated annotations at the region of interest
            grcuratedanno <- getCuratedAnno(theobject) # nolint
            gr <- .defineIntervals(grcuratedanno, regionlabel)

            if (is.na(peakpath) && isTRUE(all.equal(regionlabel, "TES"))) {

                setTermination(theobject) <- gr # nolint
                message("\t Number of ", regionlabel, ": ", length(gr))
                return(theobject)
            }else if ((isTRUE(all.equal(regionlabel, "TSS")) || #nolint
                    isTRUE(all.equal(regionlabel, "GBTES"))) &&
                !is.na(peakpath)) {

                ## Create GR for peaks
                grpeaks <- buildGR(peakpath) # nolint

                ## Overlap peaks with TSS or GBTES
                resultoverlap <- GenomicRanges::findOverlaps(grpeaks,
                        gr, ignore.strand = TRUE)

                ## Retrieve TSS coordinates having a peak
                idx <- unique(S4Vectors::subjectHits(resultoverlap))
                grresult <- gr[idx, ]

                message("\t Number of ", regionlabel, ": ", length(grresult))

                if (is.na(peaklabel) || (
                            !isTRUE(all.equal(peaklabel, "H3K27ac")) &&
                            !isTRUE(all.equal(peaklabel, "Ser5P")) &&
                            !isTRUE(all.equal(peaklabel, "Ser2P"))
                            ))
                    stop("A peak label should have been defined or the label ",
                            "provided is not valid.")
                else if (isTRUE(all.equal(regionlabel, "TSS")) &&
                        isTRUE(all.equal(peaklabel, "H3K27ac")))
                    setActivePromoters(theobject) <- grresult # nolint
                else if (isTRUE(all.equal(regionlabel, "TSS")) &&
                        isTRUE(all.equal(peaklabel, "Ser5P")))
                    setInitiationSer5PProm(theobject) <- grresult # nolint
                else if (isTRUE(all.equal(regionlabel, "GBTES")) &&
                        isTRUE(all.equal(peaklabel, "Ser2P")))
                    setElongationSer2P(theobject) <- grresult # nolint
                else
                    stop("Something went wrong in the combination of params.")

                return(theobject)
            }else
                stop("The parameters given are not correct")
        })




.overlapAndMergedTwoSamples <- function(peaks1path, peaks2pathorgr) { #nolint

    ## Create GR for each sample peaks
    grpeaks1 <- filterChromAndStrand(buildGR(peaks1path)) # nolint
    if (is.character(peaks2pathorgr))
        grpeaks2 <- filterChromAndStrand(buildGR(peaks2pathorgr)) # nolint
    else
        grpeaks2 <- peaks2pathorgr

    ## Overlap peaks of K4me3 and K27me3
    suppressWarnings(resultoverlap <- GenomicRanges::findOverlaps(grpeaks1,
      grpeaks2, ignore.strand = TRUE))

    ## Retrieving overlapping peaks
    grpeaks1 <- grpeaks1[unique(S4Vectors::queryHits(resultoverlap)), ]
    grpeaks2 <- grpeaks2[unique(S4Vectors::subjectHits(resultoverlap)), ]

    ## Intersect data to have overlapping intervals
    suppressWarnings(grmerged <- GenomicRanges::intersect(grpeaks1, grpeaks2))

    return(grmerged)
}


#' Identify Bivalent Promoters
#'
#' @description
#' This method identifies bivalent promoters by creating a merged GRanges object
#' from peaks of H3K4me3 and H3K27me3 histone marks at TSS-/+1Kb.
#'
#' @usage
#' bivalentPromoters(theobject, peakspathvec)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param peakspathvec A named vector specifying paths to the peak files of
#' H3K4me3 and H3K27me3.
#'
#' @return Returns the modified `genomicCompartments` object with bivalent
#' promoters assigned.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Identify bivalent promoters
#' bivalentPromoters(gc_obj, peakspathvec = c(H3K4me3 = "path/to/H3K4me3.gff",
#' H3K27me3 = "path/to/H3K27me3.gff"))
#' }
#'
setMethod(

        f = "bivalentPromoters",

        signature = "genomicCompartments",

        definition = function(theobject, peakspathvec) {

            ## Check the Object
            validObject(theobject)

            ## Create merged GR for peaks H3K4me3 and H3K27me3
            grmerged <- .overlapAndMergedTwoSamples(peakspathvec["H3K4me3"],
                    peakspathvec["H3K27me3"])

            message("Number of bivalent promoters: ", length(grmerged))
            ## Modify and return the object
            setBivalentProm(theobject) <- grmerged # nolint
            return(theobject)
        })


#' Identify Enhancers
#'
#' @description
#' This method identifies active or poisedd enhancers based on the overlap of
#' histone marks and ATAC-seq peaks.
#'
#' @usage
#' enhancer(theobject, peakspath1, peakspath2, peakspath3, label1, label2,
#' label3, statelabel)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param peakspath1 Path to the peaks file of the first histone mark.
#' @param peakspath2 Path to the peaks file of the second histone mark.
#' @param peakspath3 Path to the peaks file of the third histone mark or
#' ATAC-seq.
#' @param label1 Label for the first histone mark.
#' @param label2 Label for the second histone mark.
#' @param label3 Label for the third histone mark or ATAC-seq.
#' @param statelabel Label indicating the state of enhancers ('active' or
#' 'poised').
#'
#' @details
#' Active enhancers are defined as the overlap of H3K27ac/H3K4me1/ATAC-seq.
#' Poised enhancers are defined as the overlap of H3K27me3/H3K4me1/PRC2.
#' They should not overlap with the combination of UCSC refGene, NCBI RefSeq,
#' and GENCODE VM25.
#'
#' @return Returns the modified `genomicCompartments` object with enhancers
#' assigned.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Identify active enhancers
#' enhancer(gc_obj, peakspath1 = "path/to/H3K27ac.gff",
#'  peakspath2 = "path/to/H3K4me1.gff",
#'  peakspath3 = "path/to/ATAC-seq.gff",
#' label1 = "H3K27ac", label2 = "H3K4me1", label3 = "ATAC-seq",
#' statelabel = "active")
#' }
#'
setMethod(

        f = "enhancer",

        signature = "genomicCompartments",

        definition = function(theobject, peakspath1, peakspath2, peakspath3,
                label1, label2, label3, statelabel) {

            ## Check the Object
            validObject(theobject)

            if (!isTRUE(all.equal(statelabel, "active")) &&
                    !isTRUE(all.equal(statelabel, "poised")))
                stop("The state label should be active or poised")

            ## Define overlap between H3K27ac, H3K4me1, and ATAC-seq for active
            ## and H3K27me3/H3K4me1/PRC2 for poised
            ## H3K27ac over H3K4me1 or H3K27me3 over H3K4me1
            grfirst <- .overlapAndMergedTwoSamples(peakspath1, peakspath2)
            message(label1, "/", label2, ": ", length(grfirst))

            ## H3K27ac-H3K4me1 over ATAC-Seq for active or
            ## H3K27me3-H3K4me1 over Suz12 for poised
            grall <- .overlapAndMergedTwoSamples(peakspath3, grfirst)
            message(label1, "-", label2, "/", label3, ": ", length(grall))

            ## Determine peaks outside annotations (UCSC refGene, NCBI RefSeq,
            ## and GENCODE VM25)
            grlistgenesanno <- getGenesFullAnno(theobject) # nolint
            resultoverlap <- GenomicRanges::findOverlaps(grall,
                    grlistgenesanno, ignore.strand = TRUE)

            ## Retrieve indexes of peaks that do NOT overlap genes (enhancers)
            idx <- unique(S4Vectors::queryHits(resultoverlap))
            idxnot <- match(seq_len(length(grall)), idx)
            idxnot <- which(is.na(idxnot))
            grenh <- grall[idxnot, ]

            message("Number of ", statelabel, " enhancers: ", length(grenh))

            ## Updating corresponding slots
            if (isTRUE(all.equal(statelabel, "active")))
                setActiveEnhancer(theobject) <- grenh # nolint
            else
                setPoisedEnhancer(theobject) <- grenh # nolint

            return(theobject)
            })


#' Identify Polycomb Domains
#'
#' @description
#' This method identifies Polycomb domains based on the overlap of two
#' marks.
#'
#' @usage
#' polycombsDomains(theobject, peakspath1, peakspath2, label1, label2,
#' domainname)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param peakspath1 Path to the peaks file of the first polycomb.
#' @param peakspath2 Path to the peaks file of the second polycomb.
#' @param label1 Label for the first polycomb.
#' @param label2 Label for the second polycomb.
#' @param domainname Name or label for the identified Polycomb domains.
#'
#' @return Returns the modified `genomicCompartments` object with Polycomb
#' domains assigned.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Identify Polycomb domains
#' polycombsDomains(gc_obj, peakspath1 = "path/to/SUZ12.gff",
#'  peakspath2 = "path/to/RING1B.gff",
#'  label1 = "SUZ12", label2 = "RING1B", domainname = "PolycombDomains")
#' }
#'
setMethod(

        f = "polycombsDomains",

        signature = "genomicCompartments",

        definition = function(theobject, peakspath1, peakspath2, label1,
                label2, domainname) {

            ## Check the Object
            validObject(theobject)

            ## Define overlap between the two samples
            grall <- .overlapAndMergedTwoSamples(peakspath1, peakspath2)
            message(label1, "/", label2, ": ", length(grall))

            ## Updating corresponding slots
            setPolycombsDomain(theobject) <- grall # nolint

            return(theobject)
        })
