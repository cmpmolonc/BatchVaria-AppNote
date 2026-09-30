## 00_setup.R -- configuration, paths and helpers, sourced by every step.
## Contains no analysis.

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(S4Vectors)
    library(BatchVaria)
})

## ---- Paths -------------------------------------------------------------

## Set by run_all.R; otherwise the working directory.
if (!exists("PROJECT_ROOT")) {
    PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)
}

PATHS <- list(
    raw       = file.path(PROJECT_ROOT, "data", "raw"),
    processed = file.path(PROJECT_ROOT, "data", "processed"),
    results   = file.path(PROJECT_ROOT, "results"),
    figures   = file.path(PROJECT_ROOT, "figures"),
    logs      = file.path(PROJECT_ROOT, "logs")
)

for (p in PATHS) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

## GDC download cache (about 5 GB). Set GDC_DIR to reuse an existing one.
if (!exists("GDC_DIR")) {
    GDC_DIR <- Sys.getenv("GDC_DIR", file.path(PROJECT_ROOT, "GDCdata"))
}

## ---- Analysis configuration -------------------------------------------

CONFIG <- list(
    project = "TCGA-BRCA",

    ## Batch and biological variables (colData columns).
    batchVar    = "plate_id",
    bioVar      = "paper_BRCA_Subtype_PAM50",
    sampleType  = "TP",             # primary tumour

    ## Plates with fewer samples are dropped; 02 reports how many.
    minPlateSize = 3L,

    ## Patients sequenced more than once keep one RNA aliquot: "latest"
    ## keeps the one whose barcode sorts last (the later plate).
    replicateRule = "latest",

    ## Name of the uncorrected log-CPM assay.
    rawAssay = "raw",

    ## Evaluation model for the anova engine (Type II sums of squares).
    anovaFormula = ~ plate_id + paper_BRCA_Subtype_PAM50,

    ## Pooled summarisation, so that fraction x assay total is the absolute
    ## contribution. R/09_weighting_comparison.R compares it with the
    ## package default ("feature").
    anovaWeighting = "pooled",

    seed = 1L
)

## The four corrections: each method naive (preserve = NULL) and
## subtype-protected.
CORRECTIONS <- list(
    combat_naive = list(
        method = "combat", preserve = NULL
    ),
    combat_protected = list(
        method = "combat", preserve = CONFIG$bioVar
    ),
    limma_naive = list(
        method = "limma", preserve = NULL
    ),
    limma_protected = list(
        method = "limma", preserve = CONFIG$bioVar
    )
)

ASSAYS_ALL <- c(CONFIG$rawAssay, names(CORRECTIONS))

## ---- Checkpoint files --------------------------------------------------

FILES <- list(
    seRaw      = file.path(PATHS$raw, "brca_se_raw.rds"),
    aliquot    = file.path(PATHS$raw, "brca_aliquot.rds"),
    seAnnot    = file.path(PATHS$processed, "brca_se_annotated.rds"),
    seLogCPM   = file.path(PATHS$processed, "brca_se_logcpm.rds"),
    bvRaw      = file.path(PATHS$processed, "bv_raw.rds"),
    bvCorr     = file.path(PATHS$processed, "bv_corrected.rds"),
    bvProfiled = file.path(PATHS$processed, "bv_profiled.rds")
)

## ---- Helpers -----------------------------------------------------------

## Timestamped progress message.
say <- function(...) {
    message(format(Sys.time(), "[%H:%M:%S] "), ...)
}

## Read a checkpoint, naming the step that creates it if it is missing.
requireCheckpoint <- function(path, producedBy) {
    if (!file.exists(path)) {
        stop(
            "Missing checkpoint: ", path,
            "\nRun ", producedBy, " first.",
            call. = FALSE
        )
    }

    readRDS(path)
}

## Write a checkpoint and log its size.
saveCheckpoint <- function(object, path) {
    saveRDS(object, path)

    say(
        "wrote ", basename(path), " (",
        format(file.size(path) / 1024^2, digits = 3), " MB)"
    )

    invisible(path)
}

## Write a results table to results/<name>.csv.
writeTable <- function(df, name) {
    path <- file.path(PATHS$results, paste0(name, ".csv"))
    utils::write.csv(df, path, row.names = FALSE)
    say("wrote results/", basename(path))
    invisible(path)
}
