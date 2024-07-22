################################################################################
############################### Getter methods #################################
################################################################################


setGeneric(
 name = "getRefPeaks",

        def = function(theobject) {
            standardGeneric("getRefPeaks")
        },
        signature = "theobject")


setGeneric(
 name = "getGenesFullAnno",

        def = function(theobject) {
            standardGeneric("getGenesFullAnno")
        },
        signature = "theobject")

setGeneric(
 name = "getCuratedAnno",

        def = function(theobject) {
            standardGeneric("getCuratedAnno")
        },
        signature = "theobject")

setGeneric(
 name = "getActiveProm",

        def = function(theobject) {
            standardGeneric("getActiveProm")
        },
        signature = "theobject")

setGeneric(
 name = "getActiveEnh",

        def = function(theobject) {
            standardGeneric("getActiveEnh")
        },
        signature = "theobject")

setGeneric(
 name = "getPoisedEnh",

        def = function(theobject) {
            standardGeneric("getPoisedEnh")
        },
        signature = "theobject")

setGeneric(
 name = "getPcGDomain",

        def = function(theobject) {
            standardGeneric("getPcGDomain")
        },
        signature = "theobject")

setGeneric(
 name = "getHeteroChrom",

        def = function(theobject) {
            standardGeneric("getHeteroChrom")
        },
        signature = "theobject")

setGeneric(
 name = "getBivalentProm",

        def = function(theobject) {
            standardGeneric("getBivalentProm")
        },
        signature = "theobject")

setGeneric(
 name = "getInitiation",

        def = function(theobject) {
            standardGeneric("getInitiation")
        },
        signature = "theobject")

setGeneric(
 name = "getElongation",

        def = function(theobject) {
            standardGeneric("getElongation")
        },
        signature = "theobject")

setGeneric(
 name = "getTermination",

        def = function(theobject) {
            standardGeneric("getTermination")
        },
        signature = "theobject")

setGeneric(
 name = "getSINE",

        def = function(theobject) {
            standardGeneric("getSINE")
        },
        signature = "theobject")

setGeneric(
 name = "getLINE",

        def = function(theobject) {
            standardGeneric("getLINE")
        },
        signature = "theobject")

setGeneric(
 name = "getLTR",

        def = function(theobject) {
            standardGeneric("getLTR")
        },
        signature = "theobject")


setGeneric(
 name = "getRegionsMatrix",

        def = function(theobject) {
            standardGeneric("getRegionsMatrix")
        },
        signature = "theobject")


setGeneric(
 name = "getGlcPeakVal",

        def = function(theobject) {
            standardGeneric("getGlcPeakVal")
        },
        signature = "theobject")

setGeneric(
 name = "getActivePromWithPeaks",

        def = function(theobject) {
            standardGeneric("getActivePromWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getActiveEnhWithPeaks",

        def = function(theobject) {
            standardGeneric("getActiveEnhWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getPoisedEnhWithPeaks",

        def = function(theobject) {
            standardGeneric("getPoisedEnhWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getPcGDomainWithPeaks",

        def = function(theobject) {
            standardGeneric("getPcGDomainWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getHeteroChromWithPeaks",

        def = function(theobject) {
            standardGeneric("getHeteroChromWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getBivalentPromWithPeaks",

        def = function(theobject) {
            standardGeneric("getBivalentPromWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getInitiationWithPeaks",

        def = function(theobject) {
            standardGeneric("getInitiationWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getElongationWithPeaks",

        def = function(theobject) {
            standardGeneric("getElongationWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "getTerminationWithPeaks",

        def = function(theobject) {
            standardGeneric("getTerminationWithPeaks")
        },
        signature = "theobject")

setGeneric(
 name = "aslist",

        def = function(theobject, includerepeats) {
            standardGeneric("aslist")
        },
        signature = "theobject")

################################################################################
############################### Setter methods #################################
################################################################################


setGeneric(
 name = "setCurratedAnno<-",

        def = function(theobject, value) {
            standardGeneric("setCurratedAnno<-")
        },
        signature = "theobject")


setGeneric(
 name = "setActivePromoters<-",

        def = function(theobject, value) {
            standardGeneric("setActivePromoters<-")
        },
        signature = "theobject")

setGeneric(
 name = "setInitiationSer5PProm<-",

        def = function(theobject, value) {
            standardGeneric("setInitiationSer5PProm<-")
        },
        signature = "theobject")


setGeneric(
 name = "setElongationSer2P<-",

        def = function(theobject, value) {
            standardGeneric("setElongationSer2P<-")
        },
        signature = "theobject")

setGeneric(
 name = "setTermination<-",

        def = function(theobject, value) {
            standardGeneric("setTermination<-")
        },
        signature = "theobject")

setGeneric(

        name = "setBivalentProm<-",

        def = function(theobject, value) {
            standardGeneric("setBivalentProm<-")
        },
        signature = "theobject")


setGeneric(

        name = "setActiveEnhancer<-",

        def = function(theobject, value) {
            standardGeneric("setActiveEnhancer<-")
        },
        signature = "theobject")


setGeneric(

        name = "setPoisedEnhancer<-",

        def = function(theobject, value) {
            standardGeneric("setPoisedEnhancer<-")
        },
        signature = "theobject")


setGeneric(

        name = "setPolycombsDomain<-",

        def = function(theobject, value) {
            standardGeneric("setPolycombsDomain<-")
        },
        signature = "theobject")


setGeneric(

        name = "setHeterochromatin<-",

        def = function(theobject, value) {
            standardGeneric("setHeterochromatin<-")
        },
        signature = "theobject")


setGeneric(

        name = "setSINE<-",

        def = function(theobject, value) {
            standardGeneric("setSINE<-")
        },
        signature = "theobject")


setGeneric(

        name = "setLINE<-",

        def = function(theobject, value) {
            standardGeneric("setLINE<-")
        },
        signature = "theobject")


setGeneric(

        name = "setLTR<-",

        def = function(theobject, value) {
            standardGeneric("setLTR<-")
        },
        signature = "theobject")


setGeneric(

        name = "setRegionsMatrix<-",

        def = function(theobject, value) {
            standardGeneric("setRegionsMatrix<-")
        },
        signature = "theobject")


setGeneric(

        name = "setGlcPeakVal<-",

        def = function(theobject, value) {
            standardGeneric("setGlcPeakVal<-")
        },
        signature = "theobject")


setGeneric(

        name = "setActivePromWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setActivePromWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setActiveEnhWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setActiveEnhWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setPoisedEnhWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setPoisedEnhWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setPcGDomainWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setPcGDomainWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setHeteroChromWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setHeteroChromWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setBivalentPromWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setBivalentPromWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setInitiationWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setInitiationWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setElongationWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setElongationWithPeaks<-")
        },
        signature = "theobject")


setGeneric(

        name = "setTerminationWithPeaks<-",

        def = function(theobject, value) {
            standardGeneric("setTerminationWithPeaks<-")
        },
        signature = "theobject")


################################################################################
############################### genomicCompartments methods ####################
################################################################################


setGeneric(

        name = "annoCuration",

        def = function(theobject, geneannovec, nbcpu = 8) {
            standardGeneric("annoCuration")
        },
        signature = "theobject")


setGeneric(

        name = "overlapOnGenes",

        def = function(theobject, regionlabel, peaklabel = NA, peakpath = NA) {
            standardGeneric("overlapOnGenes")
        },
        signature = "theobject")

setGeneric(

        name = "bivalentPromoters",

        def = function(theobject, peakspathvec) {
            standardGeneric("bivalentPromoters")
        },
        signature = "theobject")

setGeneric(

        name = "enhancer",

        def = function(theobject, peakspath1, peakspath2, peakspath3,
                label1, label2, label3, statelabel) {
            standardGeneric("enhancer")
        },
        signature = "theobject")


setGeneric(

        name = "polycombsDomains",

        def = function(theobject, peakspath1, peakspath2, label1, label2,
            domainname) {
            standardGeneric("polycombsDomains")
        },
        signature = "theobject")


setGeneric(

        name = "upsetDiagram",

        def = function(theobject, outputfolder, includerepeats,
            plotupset = FALSE) {
            standardGeneric("upsetDiagram")
        },
        signature = "theobject")


setGeneric(

        name = "matrixForUpset",

        def = function(theobject, includerepeats, glcpeakvalues) {
            standardGeneric("matrixForUpset")
        },
        signature = "theobject")


setGeneric(

        name = "retrieveGlcPeakVal",

        def = function(theobject, includerepeats, bwpath) {
            standardGeneric("retrieveGlcPeakVal")
        },
        signature = "theobject")


setGeneric(

        name = "boxplotGlcnacLevels",

        def = function(theobject, bwpath, outputfolder, includerepeats,
            plotviolin = FALSE) {
            standardGeneric("boxplotGlcnacLevels")
        },
        signature = "theobject")


setGeneric(

        name = "outputGlcPeaksCoordPerCompartment",

        def = function(theobject, outputfolder, glcpeakspath, includerepeats) {
            standardGeneric("outputGlcPeaksCoordPerCompartment")
        },
        signature = "theobject")


setGeneric(

        name = "violinplotExpression",

        def = function(theobject, outfold, includerep, rnacounts,
                refseqpath, geneidtab, countsannotype, plotviolin = FALSE) {
            standardGeneric("violinplotExpression")
        },
        signature = "theobject")


setGeneric(

        name = "extractCompCoordWithPeak",

        def = function(theobject, outfold, includerep = FALSE) {
            standardGeneric("extractCompCoordWithPeak")
        },
        signature = "theobject")


################################################################################
############################### other methods ####################
################################################################################


setGeneric(

        name = "complexUpsetDiagram",

        def = function(theobject, includerepeats, outfold) {
            standardGeneric("complexUpsetDiagram")
        },
        signature = "theobject")
