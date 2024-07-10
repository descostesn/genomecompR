setReplaceMethod(
        f = "setCurratedAnno",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@curatedAnno <- value
            validObject(theobject)
            return(theobject)
        })

setReplaceMethod(
        f = "setActivePromoters",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@activeProm <- value
            validObject(theobject)
            return(theobject)
        })

setReplaceMethod(
        f = "setInitiationSer5PProm",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@initiation <- value
            validObject(theobject)
            return(theobject)
        })

setReplaceMethod(
        f = "setElongationSer2P",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@elongation <- value
            validObject(theobject)
            return(theobject)
        })

setReplaceMethod(
        f = "setTermination",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@termination <- value
            validObject(theobject)
            return(theobject)
        })

setReplaceMethod(
        f = "setBivalentProm",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@bivalentProm <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setActiveEnhancer",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@activeEnh <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setPoisedEnhancer",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@poisedEnh <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setPolycombsDomain",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@PcGDomain <- value
            validObject(theobject)
            return(theobject)
        })

#' Allocates heterochromatic region coordinates
#'
#' @description
#' This method stores the heterochromatin region coordinates in an object of
#' class `genomicCompartments`.
#'
#' @usage
#' setHeterochromatin(obj) <- value
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param value GRanges object of heterochromatic region coordinates.
#'
#' @return The updated object of class `genomicCompartments` with the new
#' heterochromatin regions.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- new("genomicCompartments")
#' # Define new heterochromatin regions
#' new_heterochromatin <- GRanges(seqnames = "chr1",
#' ranges = IRanges(start = 1, end = 1000000))
#' # Set the heterochromatin regions
#' setHeterochromatin(gc_obj) <- new_heterochromatin
#' }
#'
setReplaceMethod(
        f = "setHeterochromatin",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@heteroChrom <- value
            validObject(theobject)
            return(theobject)
        })

#' Allocates SINE coordinates
#'
#' @description
#' This method stores the SINE coordinates in an object of class
#' `genomicCompartments`.
#'
#' @usage
#' setSINE(obj) <- value
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param value GRanges object of SINE region coordinates.
#'
#' @return The updated object of class `genomicCompartments` with the new SINE
#' regions.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- new("genomicCompartments")
#' # Define new SINE regions
#' new_SINE <- GRanges(seqnames = "chr1",
#' ranges = IRanges(start = 1, end = 1000000))
#' # Set the SINE regions
#' setSINE(gc_obj) <- new_SINE
#' }
#'
setReplaceMethod(
        f = "setSINE",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@SINE <- value
            validObject(theobject)
            return(theobject)
        })

#' Allocate LINE coordinates
#'
#' @description
#' This method stores the LINE coordinates in an object of
#' class `genomicCompartments`.
#'
#' @usage
#' setLINE(obj) <- value
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param value GRanges object of LINE region coordinates.
#'
#' @return The updated object of class `genomicCompartments` with the new LINE
#' regions.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- new("genomicCompartments")
#' # Define new LINE regions
#' new_LINE <- GRanges(seqnames = "chr1",
#' ranges = IRanges(start = 1, end = 1000000))
#' # Set the LINE regions
#' setLINE(gc_obj) <- new_LINE
#' }
#'
setReplaceMethod(
        f = "setLINE",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@LINE <- value
            validObject(theobject)
            return(theobject)
        })

#' Allocate LTR coordinates
#'
#' @description
#' This method stores the heterochromatin region coordinates in an object of
#' class `genomicCompartments`.
#'
#' @usage
#' setLTR(obj) <- value
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param value GRanges object of LTR region coordinates.
#'
#' @return The updated object of class `genomicCompartments` with the new LTR
#' regions.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- new("genomicCompartments")
#' # Define new LTR regions
#' new_LTR <- GRanges(seqnames = "chr1",
#' ranges = IRanges(start = 1, end = 1000000))
#' # Set the LTR regions
#' setLTR(gc_obj) <- new_LINE
#' }
#'
setReplaceMethod(
        f = "setLTR",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@LTR <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setRegionsMatrix",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@regionsMatrix <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setGlcPeakVal",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@glcPeakVals <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setActivePromWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@activePromWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setActiveEnhWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@activeEnhWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setPoisedEnhWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@poisedEnhWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setPcGDomainWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@PcGDomainWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setHeteroChromWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@heteroChromWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setBivalentPromWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@bivalentPromWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setInitiationWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@initiationWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setElongationWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@elongationWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })


setReplaceMethod(
        f = "setTerminationWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject, value) {
            theobject@terminationWithPeaks <- value
            validObject(theobject)
            return(theobject)
        })
