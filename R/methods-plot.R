#################################
## matrixForUpset
#################################


.matPerPeakComp <- function(regionsperpeaklist, glcpeakvalues) { # nolint

    regionslist <- mapply(function(currentpeakregions, currentcoordinates,
                    glcpeakval) {

                namescategories <- names(glcpeakval)
                isvec <- rep(FALSE, length(namescategories))
                names(isvec) <- namescategories
                currentpeakregions <- unique(currentpeakregions)
                isvec[currentpeakregions] <- TRUE

                valpeak <- unique(unlist(lapply(glcpeakval[isvec],
                                        function(x, coord) {
                                            x[which(names(x) == coord)]},
                                            currentcoordinates)))

                if (!isTRUE(all.equal(length(valpeak),1)))
                    stop("Several values for the same peak found in ",
                            ".matPerPeakComp")

                return(c(isvec, "valpeak"=valpeak))
            }, regionsperpeaklist, names(regionsperpeaklist),
            MoreArgs = list(glcpeakvalues), SIMPLIFY = FALSE)

    regionsmatrix <- do.call(rbind, regionslist)

    if (!isTRUE(all.equal(nrow(regionsmatrix), length(regionsperpeaklist))) ||
            !isTRUE(all.equal(ncol(regionsmatrix), length(glcpeakvalues) + 1)))
        stop("regionsmatrix does not have the correct dimensions.")

    return(regionsmatrix)
}

#' Create Matrix for UpSet Plot
#'
#' @description
#' This method creates a matrix for generating an UpSet plot, indicating
#' overlaps between glc peaks and genomic compartments.
#'
#' @usage
#' matrixForUpset(theobject, includerepeats, glcpeakvalues)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param includerepeats Logical indicating whether to include repeats in the
#' analysis.
#' @param includeenhancers Logical indicating whether to include enhancers in
#' the analysis.
#' @param glcpeakvalues Numeric vector of glc peak values.
#'
#' @return Returns the modified `genomicCompartments` object with the regions
#' matrix assigned.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Create matrix for UpSet plot
#' matrixForUpset(gc_obj, includerepeats = FALSE, glcpeakvalues = c(1, 2, 3))
#' }
#'
setMethod(

        f = "matrixForUpset",

        signature = "genomicCompartments",

        definition = function(theobject, includerepeats, includeenhancers,
            glcpeakvalues) {

            ## Check the Object
            validObject(theobject)

            ## Perform the overlap
            message("Performing overlap")
            resultoverlap <- .overlapGlucnacComp(theobject, includerepeats, # nolint
                includeenhancers)

            ## Creata a list with each overlapping compartments for each peak
            complist <- aslist(theobject, includerepeats, includeenhancers) # nolint
            compnamevec <- names(complist)
            subjecthitsnames <- compnamevec[
                    S4Vectors::subjectHits(resultoverlap)]
            regionsperpeaklist <- split(subjecthitsnames,
                    S4Vectors::queryHits(resultoverlap))

            ## Indicating the peak coordinates as the names of the elements of
            ## regionsperpeaklist
            querygr <- getRefPeaks(theobject) # nolint
            coordvec <- paste(GenomicRanges::seqnames(querygr),
                GenomicRanges::start(querygr), GenomicRanges::end(querygr),
                sep = "-")
            names(regionsperpeaklist) <- coordvec[
                    as.numeric(names(regionsperpeaklist))]

            ## Create a matrix with each row corresponding to a glc peak
            ## and each column to a genomic compartments. TRUE value indicates
            ## an overlap.
            regionsmatrix <- .matPerPeakComp(regionsperpeaklist, glcpeakvalues)
            setRegionsMatrix(theobject) <- regionsmatrix # nolint

            return(theobject)
        })




#################################
## upsetDiagram
#################################


.plotUpset <- function(regionsmatrix, outputfolder, plotupset) { # nolint

    ## Removing column containing the peak values and converting to logical
    idx <- which(colnames(regionsmatrix) == "valpeak")
    regionsmatrix <- regionsmatrix[, -idx]
    namesvec <- rownames(regionsmatrix)
    regionsmatrix <- apply(regionsmatrix, 2, as.logical)
    rownames(regionsmatrix) <- namesvec

    namescolregions <- colnames(regionsmatrix)
    res <- tibble::tibble(anno = lapply(seq_len(nrow(regionsmatrix)),
                    function(i) namescolregions[regionsmatrix[i, ]]))
    g <- ggplot2::ggplot(res, aes_(x = ~anno)) + geom_bar() + # nolint
            xlab(NULL) + ylab(NULL) + theme_minimal() + # nolint
            ggupset::scale_x_upset(n_intersections = 20, order_by = "freq")

    if (plotupset) {
        dev.new()
        print(g)
    } else {
        ggplot2::ggsave(filename = "upSetGenomicCompartment.pdf", plot = g,
            device = "pdf", path = outputfolder)
    }
}

#' Generate UpSet Diagram
#'
#' @description
#' This method generates an UpSet diagram based on the overlap matrix of glc
#' peaks and genomic compartments.
#'
#' @usage
#' upsetDiagram(theobject, outputfolder, includerepeats, plotupset = FALSE)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param outputfolder Character string specifying the output folder where the
#' UpSet diagram will be saved.
#' @param includerepeats Logical indicating whether to include repeats in the
#' analysis.
#' @param plotupset Logical indicating whether to plot on device the upset plot
#' (`TRUE`) or not (`FALSE`).
#'
#' @return Does not return anything but saves the plot to outpufolder.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Generate UpSet diagram
#' upsetDiagram(gc_obj, outputfolder = "results/upset_plots",
#' includerepeats = FALSE, plotupset = TRUE)
#' }
#'
setMethod(

        f = "upsetDiagram",

        signature = "genomicCompartments",

        definition = function(theobject, outputfolder, includerepeats,
            plotupset = FALSE) {

            ## Check the Object
            validObject(theobject)

            ## Building the overlap matrix
            glcpeakvalues <- getGlcPeakVal(theobject) # nolint
            theobject <- matrixForUpset(theobject, includerepeats,
                    glcpeakvalues)

            ## Generating the upset diagram
            .plotUpset(getRegionsMatrix(theobject), outputfolder, plotupset) # nolint
        })



#################################
## retrieveGlcPeakVal
#################################


.retrieveBindingVal <- function(bwpath, complist) { # nolint

    message("Retrieving binding values (this might take a while):")

    glcnacbw <- rtracklayer::BigWigFile(bwpath)

    bindingvallist <- mapply(function(currentcomp, currentname, refbw) {

                message("\t ", currentname)

                if (!isTRUE(all.equal(length(currentcomp), 0))) {
                    mat_tmp <- rtracklayer::summary(refbw,
                            which = currentcomp,
                            size = 10000,
                            type = "mean", as = "matrix")
                    rownames(mat_tmp) <- paste(seqnames(currentcomp), # nolint
                            start(currentcomp), end(currentcomp), sep = "-")
                    return(mat_tmp)
                }else {
                    return(NULL)
                }}, complist, names(complist), MoreArgs = list(glcnacbw),
            SIMPLIFY = FALSE)

    idxnull <-  which(sapply(bindingvallist, is.null))

    if (!isTRUE(all.equal(length(idxnull), 0)))
        bindingvallist <- bindingvallist[-idxnull]

    message("Computing mean levels.")
    mean_list <- lapply(bindingvallist,
            function(mat) return(apply(mat, MARGIN = 1, mean)))
    rm(bindingvallist)
    gc()
    return(mean_list)
}

.overlapByComp <- function(complist, querygr) { # nolint

    message("Keeping compartment intervals having a glcnac peak only")
    glclist <- mapply(function(currentcompgr, currentcompname, glcpeaksgr) {

                message("\t ", currentcompname)
                result <- GenomicRanges::findOverlaps(glcpeaksgr, currentcompgr,
                        ignore.strand = TRUE)
                idx <- unique(S4Vectors::queryHits(result))
                return(glcpeaksgr[idx,])
            }, complist, names(complist), MoreArgs=list(querygr), 
            SIMPLIFY=FALSE)
    return(glclist)
}


#' Retrieve Glc Peak Values
#'
#' @description
#' This method retrieves glc peak values associated with genomic
#' compartments.
#'
#' @usage
#' retrieveGlcPeakVal(theobject, includerepeats, bwpath)
#' 
#' @param theobject An object of class `genomicCompartments`.
#' @param includerepeats Logical indicating whether to include repeats in the
#' analysis.
#' @param bwpath Character string specifying the path to the bigWig file
#' containing glc peak values.
#'
#' @return Returns the modified `genomicCompartments` object with glc peak
#' values assigned.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Retrieve glc peak values
#' gc_obj <- retrieveGlcPeakVal(gc_obj, includerepeats = FALSE,
#' bwpath = "path/to/bigWig/files")
#' }
#'
setMethod(

        f = "retrieveGlcPeakVal",

        signature = "genomicCompartments",

        definition = function(theobject, includerepeats, bwpath) {

            ## Check the Object
            validObject(theobject)

            ## Considering overlapping glcnac peaks only
            ## Retrieving lists of GR for glcnac peaks and compartments
            querygr <- getRefPeaks(theobject) # nolint
            complist <- aslist(theobject, includerepeats) # nolint

            ## Retrieving glcnac peak on each compartments
            glclist <- .overlapByComp(complist, querygr)

            ## Removing compartments without overlap
            idx <-  which(lengths(glclist) == 0)
            if (!isTRUE(all.equal(length(idx), 0))) {
                message("Number of compartments without overlap: ",
                  length(idx))
                message("\t Removing: ", paste(names(idx), collapse = " "))
                glclist <-  glclist[-idx]
            }

            mean_list <- .retrieveBindingVal(bwpath, glclist)

            setGlcPeakVal(theobject) <- mean_list # nolint

            return(theobject)
        })




#################################
## violinPlotExpression
#################################


.violinPlot <- function(mean_list, outputfilename, outputfolder, plotviolin) { # nolint

    message("Building data frame for violin plot")
    namesvec <- unlist(
            mapply(
                    function(currentmeanvec, currentname) {
                        return(rep(currentname, length(currentmeanvec)))
                    }, mean_list, names(mean_list), SIMPLIFY = FALSE))
    bindingval <- unlist(mean_list)
    df <- data.frame(val = log10(bindingval), comp = namesvec)

    message("Generating violin plot")
    g <- ggplot2::ggplot(df, ggplot2::aes(val, comp, colour = factor(comp))) + # nolint
            ggplot2::geom_violin(trim = FALSE) + ggplot2::coord_flip() +
            ggplot2::theme_classic() +
            ggplot2::theme(panel.background = ggplot2::element_rect(
                fill = 'white', colour = 'black'), # nolint
                plot.title = element_text(size = 20, face = "bold", # nolint
                            margin = ggplot2::margin(10, 0, 10, 0)),
                    axis.text.x = element_text(angle = 90, size = 10))
    g <- g + ggplot2::geom_boxplot(width = 0.1, outlier.colour = "NA")
    g <- g + ggplot2::geom_jitter(ggplot2::aes(color = factor(comp)), # nolint
        alpha = 0.2)

    if (plotviolin) {
        dev.new()
        print(g)
    }else {
        ggplot2::ggsave(filename = outputfilename, plot = g, device = "pdf",
            path = outputfolder)
    }
}

.retrieveCompartmentCoord <- function(currentname, theobject) { # nolint

    coordgr <- switch(currentname,
            "activeProm" = getActivePromWithPeaks(theobject), # nolint
            "activeEnh" = getActiveEnhWithPeaks(theobject), # nolint
            "poisedEnh" = getPoisedEnhWithPeaks(theobject), # nolint
            "PcGDomain" = getPcGDomainWithPeaks(theobject), # nolint
            "heteroChrom" = getHeteroChromWithPeaks(theobject), # nolint
            "bivalentProm" = getBivalentPromWithPeaks(theobject), # nolint
            "initiation" = getInitiationWithPeaks(theobject), # nolint
            "elongation" = getElongationWithPeaks(theobject), # nolint
            "termination" = getTerminationWithPeaks(theobject)) # nolint
    return(coordgr)
}

.returnAnnoGR <- function(idxanno, currentname, annogr) { # nolint

    if (isTRUE(all.equal(length(idxanno), 0))) {
        message("No genes were found for ", currentname)
        return(NULL)
    } else {
        return(annogr[idxanno, ])
    }
}

.retrieveAnnoGR <- function(currentname, compcoordgr, curatedgr, allgr) { # nolint

    if (isTRUE(all.equal(currentname, "activeProm")) ||
        isTRUE(all.equal(currentname, "initiation")) ||
        isTRUE(all.equal(currentname, "elongation")) ||
        isTRUE(all.equal(currentname, "termination"))) {

        res <- GenomicRanges::findOverlaps(compcoordgr, curatedgr)
        idxanno <- unique(S4Vectors::subjectHits(res))
        annogr <- .returnAnnoGR(idxanno, currentname, curatedgr)

    } else if (isTRUE(all.equal(currentname, "PcGDomain")) ||
        isTRUE(all.equal(currentname, "heteroChrom")) ||
        isTRUE(all.equal(currentname, "bivalentProm"))) {

        res <- GenomicRanges::findOverlaps(compcoordgr, allgr)
        idxanno <- unique(S4Vectors::subjectHits(res))
        annogr <- .returnAnnoGR(idxanno, currentname, allgr)

    } else {
        message("\t\t closest feature")
        idxanno <- GenomicRanges::nearest(compcoordgr, allgr, 
            ignore.strand = TRUE)
        annogr <- .returnAnnoGR(idxanno, currentname, allgr)
    }

    return(annogr)
}

.checkIdx <- function(data1, data2) { # nolint

    idx <- unique(match(data1, data2))
    idxna <- which(is.na(idx))

    if (!isTRUE(all.equal(length(idxna), 0)))
        stop("Some genes coord were not retrieved in violinplotExpression.")

    return(idx)
}

.retrieveGeneCounts <- function(countsids, rnacount, refseqgenelength, # nolint
    compcoordgr) {

    idx <- match(countsids, rnacount[, 1])
    idxna <- which(is.na(idx))

    if (!isTRUE(all.equal(length(idxna), 0))) {
        message(
            "Failing to retrieve ", length(idxna), "/",
            length(idx), " gene counts "
        )
        idx <- idx[-idxna]
        refseqgenelength <- refseqgenelength[-idxna]
    }

    if (!isTRUE(all.equal(length(idx), 0))) {
        countvec <- rnacount[idx, 2]
        x <- countvec / refseqgenelength

        if (!isTRUE(all.equal(sum(x), 0)))
            countsvectpm <- (x * 1e6) / sum(x)
        else
            countsvectpm <- x
        message("\t\t Total: ", length(countsvectpm), "/",
            length(compcoordgr), " gene expressions retrieved\n\n")
        return(countsvectpm)
    } else {
        message("No counts were retrieved")
        return(NULL)
    }
}

.tpmValuesRNASeq <- function(compnamesvec, curatedannogr, allannogr, theobject, # nolint
    coordrefseq, refseq, geneidtab, outfoldgene, rnacount, countsannotype) {

        countsvectpmlist <- sapply(compnamesvec, function(currentname,
            curatedgr, allgr, theobject, coordrefseq, refseq, geneidtab,
            outfoldgene, rnacount) {

            message("\t Retrieving genes associated to ", currentname)
            compcoordgr <- unique(.retrieveCompartmentCoord(currentname,
                theobject))
            annogr <- .retrieveAnnoGR(currentname, compcoordgr, curatedgr,
                allgr)

            if (S4Vectors::isEmpty(annogr))
                return(NULL)

            message("\t Retrieved ", length(annogr), " annotations for ",
                length(compcoordgr), " compartments")

            ## Retrieving the genes coord and ID
            coord <- paste(GenomeInfoDb::seqnames(annogr), start(annogr),
                end(annogr), sep = ":")
            idx <- .checkIdx(coord, coordrefseq)

            refseqselectionid <- gsub("\\.[0-9]+", "", refseq[idx, 3],
                perl = TRUE)
            refseqgenelength <- refseq[idx, 5] - refseq[idx, 4]

            ## Writing refseq gff to output file
            write.table(refseq[idx, ], file = file.path(outfoldgene, paste0(
                    currentname, "-genes.gff")), sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)

            ## Matching the gene ID to entrezID
            idx <- match(refseqselectionid, geneidtab[, 1])
            idxna <- which(is.na(idx))

            if (!isTRUE(all.equal(length(idxna), 0))) {
                if (isTRUE(all.equal(countsannotype, "entrez")))
                    countsids <- geneidtab[idx[-idxna], 2]
                else
                   countsids <- geneidtab[idx[-idxna], 3]

                refseqgenelength <- refseqgenelength[-idxna]
            } else {
                if (isTRUE(all.equal(countsannotype, "entrez")))
                    countsids <- geneidtab[idx, 2]
                else
                   countsids <- geneidtab[idx, 3]
            }

            ## Retrieving the gene counts
            return(.retrieveGeneCounts(countsids, rnacount,
                refseqgenelength, compcoordgr))

        }, curatedannogr, allannogr, theobject, coordrefseq,
        refseq, geneidtab, outfoldgene, rnacount
    )
    return(countsvectpmlist)
}

#' Generate Violin Plot for Gene Expression
#'
#' @description
#' This method generates violin plots for gene expression across genomic
#' compartments.
#'
#' @usage
#' violinplotExpression(theobject, outfold, includerep, rnacounts, refseqpath,
#'  geneidtab, countsannotype, plotviolin = FALSE)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param outfold Character string specifying the output folder path for saving
#' results.
#' @param includerep Logical indicating whether to include repeats in the
#' analysis.
#' @param rnacounts Character string specifying the path to the RNA-seq count
#' table.
#' @param refseqpath Character string specifying the path to the original
#' RefSeq annotation file.
#' @param geneidtab Table containing genes ID.
#' @param countsannotype Character string specifying the type of counts (e.g.,
#' entrez or not) used for expression values.
#' @param plotviolin Logical indicating whether to plot on current device
#' (default is FALSE).
#'
#' @return Save the violin plot to outfold.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Generate violin plots for gene expression
#' violinplotExpression(gc_obj, outfold = "path/to/output_folder",
#'  includerep = TRUE, rnacounts = "path/to/rna_counts",
#'  refseqpath = "path/to/refseq_annotation",
#'  geneidtab = geneidtab, countsannotype = "entrez", plotviolin = TRUE)
#' }
#'
setMethod(
    f = "violinplotExpression",

    signature = "genomicCompartments",

    definition = function(theobject, outfold, includerep, rnacounts,
                          refseqpath, geneidtab, countsannotype,
                          plotviolin = FALSE) {

        ## Check the Object
        validObject(theobject)

        ## Retrieve refseq GR for curated and all refseq anno
        curatedannogr <- getCuratedAnno(theobject) # nolint
        allannogr <- getGenesFullAnno(theobject)[[3]] # nolint

        ## Read refseq orginal file and build coord to use for matching
        refseq <- read.delim(refseqpath, header = FALSE)
        coordrefseq <- paste(refseq$V1, refseq$V4, refseq$V5, sep = ":")

        ## Overlaping each compartment with a peak with refseq GR
        ## Using the refseq GFF to retrieve gene names
        compnamesvec <- names(getGlcPeakVal(theobject)) # nolint

        ## Reading the count table
        rnacount <- read.delim(rnacounts, header = TRUE)

        ## Defining output folder for writing the genes GFF
        outfoldgene <- file.path(outfold, "genes_overlap_compWithPeaks_GFF")
        if (!file.exists(outfoldgene)) {
            dir.create(outfoldgene, recursive = TRUE)
        }

        ## Retrieve TPM values for expression
        countsvectpmlist <- .tpmValuesRNASeq(compnamesvec, curatedannogr,
            allannogr, theobject, coordrefseq, refseq, geneidtab, outfoldgene,
            rnacount, countsannotype)

        ## Remove components without counts
        idxnull <- which(sapply(countsvectpmlist,is.null))
        if (!isTRUE(all.equal(length(idxnull), 0)))
            countsvectpmlist <- countsvectpmlist[-idxnull]

        ## Generating the violin plot
        .violinPlot(countsvectpmlist, "ExpressionGenes", outfold, plotviolin)
    }
)




#################################
## boxplotGlcnacLevels
#################################

#' Generate Boxplot of Glcnac Levels
#'
#' @description
#' This method generates boxplots of Glcnac levels across the different genomic
#' compartments stored in the given genomicCompartments object.
#'
#' @usage
#' boxplotGlcnacLevels(theobject, bwpath, outputfolder, includerepeats,
#'            plotviolin = FALSE)
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param bwpath Character string specifying the path to the BigWig file
#' containing Glcnac levels.
#' @param outputfolder Character string specifying the output folder path for
#' saving results.
#' @param includerepeats Logical indicating whether to include repeats in the
#' analysis.
#' @param plotviolin Logical indicating whether to plot to current device
#' (default is FALSE).
#'
#' @return Save the boxplot to outputfolder.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Generate boxplots of Glcnac levels
#' gc_obj <- boxplotGlcnacLevels(gc_obj, bwpath = "path/to/glcnac_levels.bw",
#'                               outputfolder = "path/to/output_folder",
#'                               includerepeats = FALSE, plotviolin = FALSE)
#' }
#'
setMethod(

        f = "boxplotGlcnacLevels",

        signature = "genomicCompartments",

        definition = function(theobject, bwpath, outputfolder, includerepeats,
            plotviolin = FALSE) {

            ## Check the Object
            validObject(theobject)

            ## PART 1: Considering all annotations
            ## Retrieve glc values from bw for each compartment
            complist <- aslist(theobject, includerepeats)
            mean_list <- .retrieveBindingVal(bwpath, complist)

            ## Generating boxplot
            .violinPlot(mean_list, "violin-allcompartments-meanglc.pdf",
                outputfolder, plotviolin)

            ## PART 2: Considering overlapping glcnac peaks only
            ## Retrieving lists of GR for glcnac peaks and compartments
            theobject <- retrieveGlcPeakVal(theobject, includerepeats, bwpath)
            mean_list <- getGlcPeakVal(theobject) # nolint

            ## Generating boxplot
            .violinPlot(mean_list, "violin-glcnacPeaksCompartments-meanglc.pdf",
                outputfolder, plotviolin)
        }
)




#################################
## complexUpsetDiagram
#################################

.createMatUpset <- function(theobject, glclist, includerepeats) { # nolint

    gcmatlist <- mapply(function(currentcompartment, currentglcval,
                    includerepeats) {
                return(matrixForUpset(currentcompartment, includerepeats,
                                currentglcval))
            }, theobject, glclist, MoreArgs = list(includerepeats),
            SIMPLIFY = TRUE)
    return(gcmatlist)
}

.addColName <- function(gcmatlist) { # nolint
    gcmatlist <- mapply(function(currentobject, currentname) {

                mat <- getRegionsMatrix(currentobject) # nolint
                expmat <- as.matrix(rep(currentname, nrow(mat)))
                colnames(expmat) <- "expname"
                setRegionsMatrix(currentobject) <- cbind(expmat, mat) # nolint

                return(currentobject)
            }, gcmatlist, names(gcmatlist))
    return(gcmatlist)
}

.plotComplexUpset <- function(df, components, outfold) { # nolint

    expnamevec <- unique(df$expname)
    if (!isTRUE(all.equal(length(expnamevec), 2)))
        stop("The upset plot can only be generated for two replicates")

    g <- ComplexUpset::upset(df, components,
        base_annotations = list("Glc levels" = (
            ggplot2::ggplot(mapping = ggplot2::aes(y = as.numeric(valpeak))) + # nolint
            ggplot2::geom_jitter(ggplot2::aes(color = log10(
                as.numeric(valpeak))), na.rm = TRUE,
                alpha = 0.2) +
                ggplot2::geom_violin(alpha = 0.5, na.rm = TRUE) +
                ggplot2::coord_cartesian(ylim = quantile(as.numeric(df$valpeak),
                c(0, 0.99)))),
                "Intersection size" = ComplexUpset::intersection_size(
                    counts = FALSE, mapping = aes(fill = expname)) + # nolint
                    ggplot2::scale_fill_manual(values = c("cadetblue4",
                        "chocolate3"))), width_ratio = 0.1,
                        min_size = 20)


    ggplot2::ggsave("complexUpset.pdf", plot = g, device = "pdf",
        path = outfold)
}

#' Generate Complex Upset Diagram
#'
#' @description
#' This method generates a plot showing the upset diagram of the overlaps of
#' GlcNAc peaks with different genomic compartments and the corresponding
#' average binding values for each peak.
#'
#' @usage
#' complexUpsetDiagram(theobject, includerepeats, outfold)
#' 
#' @param theobject A list of objects of class `genomicCompartments`.
#' @param outfold Character string specifying the output folder path for saving
#' the upset diagram and boxplot.
#' @param includerepeats Logical indicating whether to include repeats in the
#' analysis. Default is FALSE.
#'
#' @return No return value. The function generates and saves a complex upset
#' plot in the specified output folder.
#'
#' @examples
#' \dontrun{
#' # Create a list of genomicCompartments objects
#' gc_list <- list(gc_obj1, gc_obj2, gc_obj3)
#' # Generate complex upset diagram
#' complexUpsetDiagram(gc_list, includerepeats = FALSE,
#'  outfold = "path/to/output_folder")
#' }
#'
setMethod(

        f = "complexUpsetDiagram",

        signature = "list",

        definition = function(theobject, outfold, includerepeats = FALSE) {

            ## Retrieve the mean levels associated to each glc peak
            glclist <- lapply(theobject, getGlcPeakVal) # nolint

            ## Build the matrices in the objects to generate the upset diagrams
            ## with chipseq levels
            gcmatlist <- .createMatUpset(theobject, glclist, includerepeats)

            ## Add an expname column to each matrix
            gcmatlist <- .addColName(gcmatlist)
            mat <- getRegionsMatrix(gcmatlist[[1]]) # nolint

            ## Verify that each replicate have the same compartments
            compnameref <- colnames(mat)
            invisible(sapply(gcmatlist[-1], function(currentobject, compref) {
                                tmpmat <- getRegionsMatrix(currentobject) # nolint
                                currentcomp <- colnames(tmpmat)
                                if (!isTRUE(all.equal(currentcomp, compref)))
                                    stop("Replicates have different comps")
                                }, compnameref))

            ## Extract and combine the matrices
            invisible(sapply(gcmatlist[-1], function(currentobject) {
                                mat <<- rbind(mat,
                                getRegionsMatrix(currentobject))})) # nolint

            ## Convert the matrix to a data.frame
            df <- as.data.frame(mat)
            components <- colnames(df)[c(-1, -ncol(df))]
            df[components] <- df[components] == "1"

            ## Generating the upset plot
            ## (See https://krassowski.github.io/complex-upset/articles/
            ## Examples_R.html)
            .plotComplexUpset(df, components, outfold)
        }
)
