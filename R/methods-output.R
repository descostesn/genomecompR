.writeCoord <- function(outfold, subfoldname, resultcoord, currentname, # nolint
        extensionname) {

    outfold <- file.path(outfold, subfoldname)
    if (!file.exists(outfold)) dir.create(outfold)
    write.table(resultcoord,
            file = file.path(outfold, paste0(currentname, extensionname)),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


#' Output Glc Peaks Coordinates per Compartment
#'
#' @description
#' This method calculates the overlap of glc peaks with each compartment
#' defined in the `genomicCompartments` object.
#'
#' @usage
#' outputGlcPeaksCoordPerCompartment(theobject, outputfolder, glcpeakspath,
#'                        includerepeats)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param outputfolder Path to the folder where output files will be written.
#' @param glcpeakspath Path to the GFF file containing glc peaks
#' coordinates.
#' @param includerepeats Logical indicating whether to include repeat regions
#' in the analysis.
#'
#' @return Returns nothing explicitly; Write GFF and BED files of the
#' coordinates of the glc peaks overlapping with a specific compartment.
#'
#' @examples
#' \dontrun{
#' # Parameters
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
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#'
#' # Output glc peaks coordinates per compartment
#' outputGlcPeaksCoordPerCompartment(gc_obj, outputfolder, peakspathquery,
#' includerepeats = FALSE)
#' }
#'
setMethod(

        f = "outputGlcPeaksCoordPerCompartment",

        signature = "genomicCompartments",

        definition = function(theobject, outputfolder, glcpeakspath,
                        includerepeats) {

            ## Check the Object
            validObject(theobject)

            ## Perform the overlap of glc peaks with each compartment
            resultoverlap <- .overlapGlucnacComp(theobject, includerepeats) # nolint

            ## Retrieve the indexes of peaks per compartment
            compnamevec <- names(aslist(theobject, includerepeats)) # nolint
            subjecthitsnames <- compnamevec[S4Vectors::subjectHits(
                resultoverlap)]
            glcpeakspercomplist <- split(S4Vectors::queryHits(resultoverlap),
                subjecthitsnames)

            ## Reading the GFF of glc peaks for retrieving coordinates
            glcgff <- read.table(glcpeakspath)

            ## Writing each coordinate
            invisible(mapply(function(currentidx, currentname, wholegfftab,
                                    outfold) {

                                message("Writing ", currentname)

                                resultgff <- wholegfftab[currentidx, ]
                                .writeCoord(outfold, "glc_coord_GFF", resultgff,
                                        currentname, "-glcPeaks.gff")

                                resultbed3col <- data.frame(resultgff$V1,
                                        resultgff$V4, resultgff$V5)
                                .writeCoord(outfold, "glc_coord_Bed_chipatlas",
                                        resultbed3col, currentname,
                                        "-glcPeaks-chipatlas.bed")

                            }, glcpeakspercomplist, names(glcpeakspercomplist),
                            MoreArgs = list(glcgff, outputfolder)))
        })



#################################
## Extract compartments coordinates that have  a glc peak
#################################


.setCompWithPeaks <- function(theobject, currentcompwithpeak, compname) { # nolint

    switch(compname,
            "activeProm" =
                setActivePromWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "activeEnh" =
                setActiveEnhWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "poisedEnh" =
                setPoisedEnhWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "PcGDomain" =
                setPcGDomainWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "heteroChrom" =
                setHeteroChromWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "bivalentProm" =
                setBivalentPromWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "initiation" =
                setInitiationWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "elongation" =
                setElongationWithPeaks(theobject) <- currentcompwithpeak, # nolint
            "termination" =
                setTerminationWithPeaks(theobject) <- currentcompwithpeak) # nolint

    return(theobject)
}

.defineDF <- function(currentcompwithpeak, compname) { # nolint
    return(data.frame(seqname = seqnames(currentcompwithpeak), # nolint
                    source ="compCoordWithPeaks",
                    feature = compname,
                    start = start(currentcompwithpeak),
                    end = end(currentcompwithpeak),
                    score = 0,
                    strand = strand(currentcompwithpeak), # nolint
                    frame = ".",
                    group = "."))
}

#' Extract Compartment Coordinates with Peaks
#'
#' @usage
#' extractCompCoordWithPeak(theobject, outfold, includerep)
#'
#' @description
#' This method performs overlap analysis between reference peaks and each
#' compartment defined in the `genomicCompartments` object. It extracts
#' compartments containing peaks and outputs their coordinates in GFF format.
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param outfold Path to the folder where output files will be written.
#' @param includerep Logical indicating whether to include repeat regions in
#' the analysis. Default is FALSE.
#' @param includeenhancers Logical indicating whether to include enhancers in
#' the analysis. Default is FALSE.
#'
#' @return Returns the modified `genomicCompartments` object with compartments
#' containing peaks.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Extract compartment coordinates with peaks
#' extractCompCoordWithPeak(gc_obj, outfold = "output_folder")
#' }
#'
setMethod(

        f = "extractCompCoordWithPeak",

        signature = "genomicCompartments",

        definition = function(theobject, outfold, includerep = FALSE,
            includeenhancers = FALSE) {

            ## Check the Object
            validObject(theobject)

            ## Overlap with the different compartments and retrieve
            ## peaks per compartment type
            message("Performing overlap with each compartment")
            complist <- aslist(theobject, includerep, includeenhancers) # nolint
            querygr <- getRefPeaks(theobject) # nolint

            invisible(mapply(function(currentcomp, compname, peaksgr, outfold) {

                        message("\t\t ", compname)
                        resultoverlap <- GenomicRanges::findOverlaps(peaksgr,
                                currentcomp, ignore.strand = TRUE)
                        idxcomp <- subjectHits(resultoverlap) # nolint
                        message("\t\t", length(unique(idxcomp)), "/",
                                length(currentcomp),
                                " compartments contain a peak")

                        if (!isTRUE(all.equal(length(idxcomp), 0))) {

                            currentcompwithpeak <- currentcomp[idxcomp, ]
                            gffdf <- .defineDF(currentcompwithpeak, compname)
                            .writeCoord(outfold, "comp_coord_withPeak_GFF",
                                    gffdf, compname, "-withPeaks.gff")

                            theobject <<- .setCompWithPeaks(theobject,
                                    currentcompwithpeak, compname)
                        }else {
                                message("## No peaks were found in ", compname)
                        }
                    }, complist, names(complist), MoreArgs = list(querygr,
                            outfold)))

            return(theobject)
        })
