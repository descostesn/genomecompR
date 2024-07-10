#' Curate Gene Annotations
#'
#' @usage
#' annoCuration(theobject, geneannovec, nbcpu = 8)
#'
#' @description
#' This method curates gene annotations by performing several filtering steps
#' on the provided annotation data. It filters out non-canonical chromosomes,
#' retains known coding and non-coding genes, keeps genes with a length greater
#' than 2kb, removes duplicate annotations with the same coordinates, extends
#' the transcription start site (TSS) by 1kb upstream, and removes overlapping
#' annotations. The curated annotations are then set as the curated annotation
#' data of the `genomicCompartments` object.
#'
#' @param theobject An object of class `genomicCompartments`.
#' @param geneannovec A character vector containing the paths to gene
#' annotation files.
#' @param nbcpu Number of CPU cores to use for parallel processing
#' (default is 8).
#'
#' @return Returns the updated `genomicCompartments` object with curated gene
#' annotations.
#'
#' @examples
#' \dontrun{
#' # Create a genomicCompartments object
#' gc_obj <- genomeCompart(peakspathvec, geneannovec)
#' # Perform annotation curation
#' gc_obj <- annoCuration(gc_obj, geneannovec)
#' }
#'
setMethod(

        f = "annoCuration",

        signature = "genomicCompartments",

        definition = function(theobject, geneannovec, nbcpu = 8){

            ## Check the Object
            validObject(theobject)

            anno <- read.table(geneannovec["refseq"])
            message("Initial number of annotations: ", nrow(anno))

            ## Filtering non-canonical chromosomes
            anno <- .filterChrom(anno) # nolint
            message("Number annotations after filtering chroms: ", nrow(anno))

            ## Retrieve known coding and non-coding genes
            idx <- c(grep("NM",anno$V3), grep("NR",anno$V3))
            anno <- anno[idx,]
            message("Number coding and non-coding: ", nrow(anno))

            ## Keeping genes with length > 2kb
            idx <- which((anno$V5 - anno$V4) > 2000)
            anno <- anno[idx,]
            message("Number anno > 2kb: ", nrow(anno))

            ## Remove annotations having exactly the same coordinates keeping
            ## one copy
            annosemi <- paste(anno$V1, anno$V4, anno$V5, sep = ":")
            anno <- anno[!duplicated(annosemi), ]
            message("Number unique annotations: ", nrow(anno))

            ## Extend TSS by 1kb upstream and remove overlapping annotations
            annoextendedup <- anno
            annoextendedup$V4 <- anno$V4 - 1000

            ## Converting annotations to a list to use parallel processing
            annoextendedup_list <- apply(annoextendedup, MARGIN=1,as.list)
            annoextendedup_list <- lapply(annoextendedup_list, as.character)
            to_keep <- list()
            message("Removing overlapping annotations")
            to_keep <- parallel::mclapply(annoextendedup_list,
                    function(x, annotations) {

                        chromosome <- x[1]
                        currentbegin <- as.numeric(x[4])
                        currentend <- as.numeric(x[5])

                        overlapind_begin_inside	<- (
                                    (currentbegin < as.numeric(annotations[,
                                                                        4]) &
                                        currentend > as.numeric(annotations[,
                                                                        4])) &
                                    chromosome == annotations[, 1])
                        overlapind_end_inside	<- (
                                    (currentbegin < as.numeric(annotations[,
                                                                        5]) &
                                        currentend > as.numeric(annotations[,
                                                                        5])) &
                                    chromosome == annotations[, 1])
                        overlap_included <- (
                                    (currentbegin > as.numeric(annotations[,
                                                                        4]) &
                                        currentend < as.numeric(annotations[,
                                                                        5])) &
                                    chromosome == annotations[, 1])
                        overlapind <- (overlapind_begin_inside |
                                    overlapind_end_inside | overlap_included)

                        if (sum(overlapind) == 0)
                            return(TRUE)
                        else
                            return(FALSE)
                    }, annoextendedup, mc.cores = nbcpu)

            to_keep <- which(unlist(to_keep))
            message("Number of annotations after removal: ", length(to_keep),
                    "\n")
            anno <- anno[to_keep, ]
            granno <- buildGRNoReading(anno) # nolint

            setCurratedAnno(theobject) <- granno # nolint

            return(theobject)
        })
