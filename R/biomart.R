
## Code taken from sharedInternals.R of the CONCLUS package

#' Attempt to Connect to Ensembl Using biomaRt
#'
#' @usage
#' tryUseMart(biomart = "ensembl", dataset, host, alternativemirror)
#'
#' @description
#' This function attempts to establish a connection to the Ensembl database
#' using #' the `biomaRt` package. It tries up to 5 times before giving up and
#' raising an error. The function supports using an alternative mirror if
#' specified.
#'
#' @param biomart Character string specifying the BioMart database to use.
#' Default is "ensembl".
#' @param dataset Character string specifying the dataset to use within the
#' BioMart database.
#' @param host Character string specifying the host URL to connect to.
#' @param alternativemirror Logical value indicating whether to use an
#' alternative mirror. Default is `FALSE`.
#'
#' @return Returns a `Mart` object if the connection is successful.
#' @examples
#' \dontrun{
#' ensembl <- tryUseMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
#' host = "https://www.ensembl.org")
#' }
#'
tryUseMart <- function(biomart = "ensembl", dataset, host, # nolint
        alternativemirror) {

    c <- 1

    repeat {
        message("# Attempt ", c, "/5 # ", "Connection to Ensembl ... ")

        if (!alternativemirror)
            ensembl <- try(biomaRt::useMart(biomart, dataset = dataset,
                host = host), silent = TRUE)
        else
            ensembl <- try(biomaRt::useEnsembl(biomart, dataset = dataset,
                host = host, mirror = "useast"), silent = TRUE)

        if (isTRUE(is(ensembl, "try-error"))) {
            c <- c + 1
            error_type <- attr(ensembl, "condition")
            message(error_type$message)

            if (c > 5)
                stop("There is a problem of connexion to Ensembl for ",
                        "now. Please retry later or set ",
                        "alternativemirror=TRUE.")
        }else {
            message("Connected with success.")
            return(ensembl)
        }
    }
}


#' Attempt to Retrieve Information from Ensembl Using biomaRt
#'
#' @usage
#' tryGetBM(attributes, ensembl, values = NULL, filters = NULL)
#'
#' @description
#' This function attempts to retrieve information about genes from Ensembl using
#' the `biomaRt` package. It tries up to 5 times before giving up and raising an
#' error. The function supports specifying attributes, values, and filters for
#' the query.
#'
#' @param attributes Character vector specifying the attributes to retrieve.
#' @param ensembl `Mart` object representing the connection to the Ensembl
#' database.
#' @param values Optional. List or vector of values for the filter.
#' @param filters Optional. Character vector specifying the filters to apply.
#'
#' @return Returns a data frame with the retrieved information if the query is
#' successful.
#'
#' @examples
#' \dontrun{
#' ensembl <- tryUseMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
#' host = "https://www.ensembl.org")
#' attributes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name",
#' "start_position", "end_position")
#' gene_info <- tryGetBM(attributes = attributes, ensembl = ensembl)
#' }
#'
tryGetBM <- function(attributes, ensembl, values = NULL, filters = NULL) {# nolint

    c <- 1

    repeat {

        message("# Attempt ", c, "/5 # ",
                "Retrieving information about genes from biomaRt ...")


        if (is.null(values) && is.null(filters))
            res <- try(biomaRt::getBM(attributes = attributes, mart = ensembl),
                silent = TRUE)
        else
            res <- try(biomaRt::getBM(attributes = attributes, mart = ensembl,
                values = values, filters = filters), silent = TRUE)

        if (isTRUE(is(res, "try-error"))) {
            c <- c + 1
            error_type <- attr(res, "condition")
            message(error_type$message)

            if (c > 5)
                stop("There is a problem of connexion to Ensembl for ",
                        "now. Please retry later.")

        }else {
            message("Information retrieved with success.")
            return(res)
        }
    }
}
