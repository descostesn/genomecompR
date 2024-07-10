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
#' @slot glcPeakVals A list of glucnac peaks.
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
#' getRefPeaks(x): Get the reference glc peaks
#' getGenesFullAnno(x): Get the full gene annotations
#' getCuratedAnno(x): Get the curated annotations
#' getActiveProm(x): Get the active promoters
#' getActiveEnh(x): Get the active enhancers
#' getPoisedEnh(x): Get the poised enhancers
#' getPcGDomain(x): Get the PcG domains
#' getHeteroChrom(x): Get the heterochromatin regions
#' getBivalentProm(x): Get the bivalent promoters
#' getInitiation(x): Get the initiation regions
#' getElongation(x): Get the elongation regions
#' getTermination(x): Get the termination regions
#' getSINE(x): Get the SINE regions
#' getLINE(x): Get the LINE regions
#' getLTR(x): Get the LTR regions
#' getRegionsMatrix(x): Get the regions matrix
#' getGlcPeakVal(x): Get the glc peak values
#' getActivePromWithPeaks(x): Get the active promoters with glc peaks
#' getActiveEnhWithPeaks(x): Get the active enhancers with glc peaks
#' getPoisedEnhWithPeaks(x): Get the poised enhancers with glc peaks
#' getPcGDomainWithPeaks(x): Get the PcG domains with glc peaks
#' getHeteroChromWithPeaks(x): Get the heterochromatin regions with glc peaks
#' getBivalentPromWithPeaks(x): Get the bivalent promoters with glc peaks
#' getInitiationWithPeaks(x): Get the initiation regions with glc peaks
#' getElongationWithPeaks(x): Get the elongation regions with glc peaks
#' getTerminationWithPeaks(x): Get the termination regions with glc peaks
#'
#' @section Subsetting:
#'
#' In the following snippets, x is a genomicCompartments object.
#'
#' setCurratedAnno(x, value): Set the curated annotations
#' setActivePromoters(x, value): Set the active promoters
#' setInitiationSer5PProm(x, value): Set the initiation Ser5P promoters
#' setElongationSer2P(x, value): Set the elongation Ser2P regions
#' setTermination(x, value): Set the termination regions
#' setBivalentProm(x, value): Set the bivalent promoters
#' setActiveEnhancer(x, value): Set the active enhancers
#' setPoisedEnhancer(x, value): Set the poised enhancers
#' setPolycombsDomain(x, value): Set the PcG domains
#' setHeterochromatin(x, value): Set the heterochromatin regions
#' setSINE(x, value): Set the SINE regions
#' setLINE(x, value): Set the LINE regions
#' setLTR(x, value): Set the LTR regions
#' setRegionsMatrix(x, value): Set the regions matrix
#' setGlcPeakVal(x, value): Set the glc peak values
#' setActivePromWithPeaks(x, value): Set the active promoters with glc peaks
#' setActiveEnhWithPeaks(x, value): Set the active enhancers with glc peaks
#' setPoisedEnhWithPeaks(x, value): Set the poised enhancers with glc peaks
#' setPcGDomainWithPeaks(x, value): Set the PcG domains with glc peaks
#' setHeteroChromWithPeaks(x, value): Set the heterochromatin regions with glc peaks
#' setBivalentPromWithPeaks(x, value): Set the bivalent promoters with glc peaks
#' setInitiationWithPeaks(x, value): Set the initiation regions with glc peaks
#' setElongationWithPeaks(x, value): Set the elongation regions with glc peaks
#' setTerminationWithPeaks(x, value): Set the termination regions with glc peaks
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
