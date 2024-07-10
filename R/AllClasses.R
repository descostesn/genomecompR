#' Genomic Compartments Class
#'
#' @description
#' The `genomicCompartments` class is designed to store  various genomic data,
#' including reference peaks, gene annotations, and different types of genomic
#' regions such as promoters, enhancers, and chromatin domains.
#'
#' @rdname genomicCompartments-class
#' @aliases genomicCompartments-class genomicCompartments
#'
#' @section Constructor:
#' new("genomicCompartments",
#'            referencePeaks = grglc,
#'            geneAnnotations = grlist, ...)
#'
#' @slot referencePeaks A `GRanges` object containing reference peaks
#' (required).
#' @slot geneAnnotations A `GRangesList` object containing gene annotations
#' (required).
#' @slot curatedAnno A `GRanges` object for curated annotations.
#' @slot activeProm A `GRanges` object for active promoters defined as K27ac in
#' TSS-/+1Kb.
#' @slot activeEnh A `GRanges` object for active enhancers defined as
#' H3K27ac/H3K4me1/ATAC-seq. They should not overlap with the combination of
#' UCSC refGene, NCBI RefSeq, and GENCODE VM25.
#' @slot poisedEnh A `GRanges` object for poised enhancers defined as
#' H3K27me3/H3K4me1/PRC2. They should not overlap with the combination of UCSC
#' refGene, NCBI RefSeq, and GENCODE VM25. Note: use a non redundant set of
#' genes i.e. no overlap, remove genes closer to 1kb from one another,
#' remove genes that are too short i.e. length < 2kb.
#' @slot PcGDomain A `GRanges` object for PcG domains defined as Suz12 and
#' RING1B overlap.
#' @slot heteroChrom A `GRanges` object for heterochromatin regions (H3K9me3)
#' @slot bivalentProm A `GRanges` object for bivalent promoters defined as
#' the overlap of H3K4me3/H3K27me3.
#' @slot initiation A `GRanges` object for initiation regions defined as
#' TSS-1/+1kb overlapping with Ser5P peaks.
#' @slot elongation A `GRanges` object for elongation regions defined as
#' TSS+1kb to TES overlapping Ser2P.
#' @slot termination A `GRanges` object for termination regions defined as
#' TES+1kb
#' @slot SINE A `GRanges` object for SINE elements.
#' @slot LINE A `GRanges` object for LINE elements.
#' @slot LTR A `GRanges` object for LTR elements.
#' @slot regionsMatrix A matrix representing various regions.
#' @slot glcPeakVals A list of glc peaks.
#' @slot activePromWithPeaks A `GRanges` object for active promoters with glc
#' peaks.
#' @slot activeEnhWithPeaks A `GRanges` object for active enhancers with glc
#' peaks.
#' @slot poisedEnhWithPeaks A `GRanges` object for poised enhancers with glc
#' peaks.
#' @slot PcGDomainWithPeaks A `GRanges` object for PcG domains with glc peaks.
#' @slot heteroChromWithPeaks A `GRanges` object for heterochromatin regions
#' with glc peaks.
#' @slot bivalentPromWithPeaks A `GRanges` object for bivalent promoters with
#' glc peaks.
#' @slot initiationWithPeaks A `GRanges` object for initiation regions with
#' glc peaks.
#' @slot elongationWithPeaks A `GRanges` object for elongation regions with glc
#' peaks.
#' @slot terminationWithPeaks A `GRanges` object for termination regions with
#' glc peaks.
#'
#' 
#' @section Accessors:
#'
#' In the following snippets, x is a genomicCompartments object.
#'
#' getRefPeaks(x): Get the reference glc peaks \cr
#' getGenesFullAnno(x): Get the full gene annotations \cr
#' getCuratedAnno(x): Get the curated annotations \cr
#' getActiveProm(x): Get the active promoters \cr
#' getActiveEnh(x): Get the active enhancers \cr
#' getPoisedEnh(x): Get the poised enhancers \cr
#' getPcGDomain(x): Get the PcG domains \cr
#' getHeteroChrom(x): Get the heterochromatin regions \cr
#' getBivalentProm(x): Get the bivalent promoters \cr
#' getInitiation(x): Get the initiation regions \cr
#' getElongation(x): Get the elongation regions \cr
#' getTermination(x): Get the termination regions \cr
#' getSINE(x): Get the SINE regions \cr
#' getLINE(x): Get the LINE regions \cr
#' getLTR(x): Get the LTR regions \cr
#' getRegionsMatrix(x): Get the regions matrix \cr
#' getGlcPeakVal(x): Get the glc peak values \cr
#' getActivePromWithPeaks(x): Get the active promoters with glc peaks \cr
#' getActiveEnhWithPeaks(x): Get the active enhancers with glc peaks \cr
#' getPoisedEnhWithPeaks(x): Get the poised enhancers with glc peaks \cr
#' getPcGDomainWithPeaks(x): Get the PcG domains with glc peaks \cr
#' getHeteroChromWithPeaks(x): Get the heterochromatin regions with glc
#' peaks \cr
#' getBivalentPromWithPeaks(x): Get the bivalent promoters with glc peaks \cr
#' getInitiationWithPeaks(x): Get the initiation regions with glc peaks \cr
#' getElongationWithPeaks(x): Get the elongation regions with glc peaks \cr
#' getTerminationWithPeaks(x): Get the termination regions with glc peaks \cr
#'
#' @section Subsetting:
#'
#' In the following snippets, x is a genomicCompartments object.
#'
#' setCurratedAnno(x, value): Set the curated annotations \cr
#' setActivePromoters(x, value): Set the active promoters \cr
#' setInitiationSer5PProm(x, value): Set the initiation Ser5P promoters \cr
#' setElongationSer2P(x, value): Set the elongation Ser2P regions \cr
#' setTermination(x, value): Set the termination regions \cr
#' setBivalentProm(x, value): Set the bivalent promoters \cr
#' setActiveEnhancer(x, value): Set the active enhancers \cr
#' setPoisedEnhancer(x, value): Set the poised enhancers \cr
#' setPolycombsDomain(x, value): Set the PcG domains \cr
#' setHeterochromatin(x, value): Set the heterochromatin regions \cr
#' setSINE(x, value): Set the SINE regions \cr
#' setLINE(x, value): Set the LINE regions \cr
#' setLTR(x, value): Set the LTR regions \cr
#' setRegionsMatrix(x, value): Set the regions matrix \cr
#' setGlcPeakVal(x, value): Set the glc peak values \cr
#' setActivePromWithPeaks(x, value): Set the active promoters with glc peaks \cr
#' setActiveEnhWithPeaks(x, value): Set the active enhancers with glc peaks \cr
#' setPoisedEnhWithPeaks(x, value): Set the poised enhancers with glc peaks \cr
#' setPcGDomainWithPeaks(x, value): Set the PcG domains with glc peaks \cr
#' setHeteroChromWithPeaks(x, value): Set the heterochromatin regions with glc
#' peaks \cr
#' setBivalentPromWithPeaks(x, value): Set the bivalent promoters with glc
#' peaks \cr
#' setInitiationWithPeaks(x, value): Set the initiation regions with glc
#' peaks \cr
#' setElongationWithPeaks(x, value): Set the elongation regions with glc
#' peaks \cr
#' setTerminationWithPeaks(x, value): Set the termination regions with glc
#' peaks \cr
#'
genomicCompartments <- setClass( # nolint
        "genomicCompartments",
        slots = c(
                referencePeaks = "GRanges",
                geneAnnotations = "GRangesList",
                curatedAnno = "GRanges",
                activeProm = "GRanges",
                activeEnh = "GRanges",
                poisedEnh = "GRanges",
                PcGDomain = "GRanges",
                heteroChrom = "GRanges",
                bivalentProm = "GRanges",
                initiation = "GRanges",
                elongation = "GRanges",
                termination = "GRanges",
                SINE = "GRanges",
                LINE = "GRanges",
                LTR = "GRanges",
                regionsMatrix = "matrix",
                glcPeakVals = "list",
                activePromWithPeaks = "GRanges",
                activeEnhWithPeaks = "GRanges",
                poisedEnhWithPeaks = "GRanges",
                PcGDomainWithPeaks = "GRanges",
                heteroChromWithPeaks = "GRanges",
                bivalentPromWithPeaks = "GRanges",
                initiationWithPeaks = "GRanges",
                elongationWithPeaks = "GRanges",
                terminationWithPeaks = "GRanges"
        ),

        prototype = list(

                curatedAnno = GenomicRanges::GRanges(),
                activeProm = GenomicRanges::GRanges(),
                activeEnh = GenomicRanges::GRanges(),
                poisedEnh = GenomicRanges::GRanges(),
                PcGDomain = GenomicRanges::GRanges(),
                heteroChrom = GenomicRanges::GRanges(),
                bivalentProm = GenomicRanges::GRanges(),
                initiation = GenomicRanges::GRanges(),
                elongation = GenomicRanges::GRanges(),
                termination = GenomicRanges::GRanges(),
                SINE = GenomicRanges::GRanges(),
                LINE = GenomicRanges::GRanges(),
                LTR = GenomicRanges::GRanges(),
                regionsMatrix = matrix(),
                glcPeakVals = list(),
                activePromWithPeaks = GenomicRanges::GRanges(),
                activeEnhWithPeaks = GenomicRanges::GRanges(),
                poisedEnhWithPeaks = GenomicRanges::GRanges(),
                PcGDomainWithPeaks = GenomicRanges::GRanges(),
                heteroChromWithPeaks = GenomicRanges::GRanges(),
                bivalentPromWithPeaks = GenomicRanges::GRanges(),
                initiationWithPeaks = GenomicRanges::GRanges(),
                elongationWithPeaks = GenomicRanges::GRanges(),
                terminationWithPeaks = GenomicRanges::GRanges()
        ),

        validity = function(object) {

            refpeaks <- getRefPeaks(object)

            if (isTRUE(all.equal(nrow(refpeaks), 0)))
                stop("The reference peaks cannot be empty.")

            genesanno <- getGenesFullAnno(object)

            if (isTRUE(all.equal(nrow(genesanno), 0)))
                stop("The genes annotations cannot be empty.")
        })
