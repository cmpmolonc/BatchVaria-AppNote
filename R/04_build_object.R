## 04_build_object.R -- establish the starting object.
##
## Checks the SummarizedExperiment from 03 before corrections are run.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

se <- requireCheckpoint(FILES$seLogCPM, "R/03_preprocess.R")

set.seed(CONFIG$seed)

bv <- se

## The starting object is valid and carries the expected assay.
stopifnot(
    is(bv, "SummarizedExperiment"),
    validObject(bv),
    CONFIG$rawAssay %in% assayNames(bv),
    identical(colnames(bv), rownames(colData(bv))),
    ## basisRetention() matches features by name.
    !is.null(rownames(bv))
)

say(
    "04: starting object: ", nrow(bv), " features x ", ncol(bv),
    " samples, assay '", CONFIG$rawAssay, "'"
)

## Batch and biological variables are present and usable.
stopifnot(
    CONFIG$batchVar %in% colnames(colData(bv)),
    CONFIG$bioVar %in% colnames(colData(bv)),
    is.factor(colData(bv)[[CONFIG$batchVar]]),
    is.factor(colData(bv)[[CONFIG$bioVar]]),
    nlevels(colData(bv)[[CONFIG$batchVar]]) > 1L,
    nlevels(colData(bv)[[CONFIG$bioVar]]) > 1L
)

say(
    "04: batch '", CONFIG$batchVar, "' has ",
    nlevels(colData(bv)[[CONFIG$batchVar]]), " levels; biological variable '",
    CONFIG$bioVar, "' has ", nlevels(colData(bv)[[CONFIG$bioVar]]), " levels"
)

## Residual degrees of freedom under the evaluation model.
designMatrix <- model.matrix(
    CONFIG$anovaFormula,
    data = as.data.frame(colData(bv))
)

say(
    "04: fixed-effect design has ", ncol(designMatrix), " columns, rank ",
    qr(designMatrix)$rank, ", residual df ",
    ncol(bv) - qr(designMatrix)$rank
)

saveCheckpoint(bv, FILES$bvRaw)
say("04: done")
