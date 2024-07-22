setMethod(
        f = "getRefPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@referencePeaks)
        })

setMethod(
        f = "getGenesFullAnno",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@geneAnnotations)
        })

setMethod(
        f = "getCuratedAnno",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@curatedAnno)
        })


setMethod(
        f = "getActiveProm",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@activeProm)
        })

setMethod(
        f = "getActiveEnh",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@activeEnh)
        })

setMethod(
        f = "getPoisedEnh",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@poisedEnh)
        })

setMethod(
        f = "getPcGDomain",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@PcGDomain)
        })

setMethod(
        f = "getHeteroChrom",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@heteroChrom)
        })

setMethod(
        f = "getBivalentProm",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@bivalentProm)
        })

setMethod(
        f = "getInitiation",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@initiation)
        })

setMethod(
        f = "getElongation",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@elongation)
        })

setMethod(
        f = "getTermination",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@termination)
        })

setMethod(
        f = "getSINE",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@SINE)
        })

setMethod(
        f = "getLINE",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@LINE)
        })

setMethod(
        f = "getLTR",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@LTR)
        })


setMethod(
        f = "getRegionsMatrix",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@regionsMatrix)
        })


setMethod(
        f = "getGlcPeakVal",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@glcPeakVals)
        })


setMethod(
        f = "aslist",
        signature = "genomicCompartments",
        definition = function(theobject, includerepeats, includeenhancers) {

            compnames <- c("activeProm", "PcGDomain", "heteroChrom",
                "bivalentProm", "initiation", "elongation", "termination")

            result <- list(theobject@activeProm, theobject@PcGDomain,
                    theobject@heteroChrom, theobject@bivalentProm,
                    theobject@initiation, theobject@elongation,
                    theobject@termination)

            if (includeenhancers) {
                compnames <- c(compnames, c("activeEnh", "poisedEnh"))
                result <- c(result, list(theobject@activeEnh,
                    theobject@poisedEnh))
            }

            if (includerepeats) {
                compnames <- c(compnames, c("SINE", "LINE", "LTR"))
                result <- c(result, list(theobject@SINE, theobject@LINE,
                    theobject@LTR))
            }

            ## Removing compartments that were not defined
            idx <-  which(lengths(result) == 0)
            if (!isTRUE(all.equal(length(idx), 0))) {
                message("Number of empty compartments: ", length(idx))
                message("Removing: ", paste(compnames[idx], collapse = " "))
                result <-  result[-idx]
                compnames <- compnames[-idx]
            }

            names(result) <- compnames
            return(result)
        })


setMethod(
        f = "getActivePromWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@activePromWithPeaks)
        })

setMethod(
        f = "getActiveEnhWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@activeEnhWithPeaks)
        })

setMethod(
        f = "getPoisedEnhWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@poisedEnhWithPeaks)
        })

setMethod(
        f = "getPcGDomainWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@PcGDomainWithPeaks)
        })

setMethod(
        f = "getHeteroChromWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@heteroChromWithPeaks)
        })

setMethod(
        f = "getBivalentPromWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@bivalentPromWithPeaks)
        })

setMethod(
        f = "getInitiationWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@initiationWithPeaks)
        })

setMethod(
        f = "getElongationWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@elongationWithPeaks)
        })

setMethod(
        f = "getTerminationWithPeaks",
        signature = "genomicCompartments",
        definition = function(theobject) {
            return(theobject@terminationWithPeaks)
        })
