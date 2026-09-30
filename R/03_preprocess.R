## 03_preprocess.R -- filtering, normalisation and log transformation.
##
## Produces the uncorrected log2-CPM matrix that every later step uses. The
## gene filter is computed on the samples fixed in 02.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

suppressPackageStartupMessages({
    library(edgeR)
    library(matrixStats)
})

se <- requireCheckpoint(FILES$seAnnot, "R/02_sample_metadata.R")

say("03: starting from ", nrow(se), " features x ", ncol(se), " samples")

## ---- Raw counts ---------------------------------------------------------

## Unstranded STAR counts.
if (!"unstranded" %in% assayNames(se)) {
    stop(
        "Expected an 'unstranded' assay. Available: ",
        paste(assayNames(se), collapse = ", ")
    )
}

counts <- assay(se, "unstranded")

stopifnot(
    is.matrix(counts),
    !is.null(dimnames(counts)),
    !anyNA(counts)
)

storage.mode(counts) <- "integer"

## ---- Filter, normalise, transform ---------------------------------------

dge <- DGEList(counts = counts)

## Grouping by subtype retains genes expressed in only one subtype.
keepExpr <- filterByExpr(dge, group = colData(se)[[CONFIG$bioVar]])
say("03: ", sum(keepExpr), " of ", nrow(dge), " features pass filterByExpr")

dge <- dge[keepExpr, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

## ---- Drop invariant features -------------------------------------------

## Drop zero-variance features so every assay uses the same feature set.
geneVar <- rowVars(logCPM)
keepVar <- is.finite(geneVar) & geneVar > 0

say("03: ", sum(!keepVar), " features dropped for zero or non-finite variance")

logCPM <- logCPM[keepVar, , drop = FALSE]

## ---- Rebuild the container ---------------------------------------------

se <- se[rownames(logCPM), ]

stopifnot(
    identical(rownames(logCPM), rownames(se)),
    identical(colnames(logCPM), colnames(se))
)

## A plain SummarizedExperiment holding only the log-CPM matrix; genomic
## coordinates are not used.
se <- SummarizedExperiment(
    assays   = setNames(SimpleList(logCPM), CONFIG$rawAssay),
    colData  = colData(se),
    rowData  = rowData(se),
    metadata = metadata(se)
)

rownames(se) <- rownames(logCPM)

stopifnot(
    validObject(se),
    !is.null(rownames(se)),
    identical(rownames(se), rownames(logCPM))
)

say("03: final matrix ", nrow(se), " features x ", ncol(se), " samples")

writeTable(
    data.frame(
        step = c("input", "filterByExpr", "variance_filter"),
        n_features = c(nrow(counts), sum(keepExpr), nrow(se))
    ),
    "03_feature_filtering"
)

saveCheckpoint(se, FILES$seLogCPM)
say("03: done")
