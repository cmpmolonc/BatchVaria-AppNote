## 01_download_tcga.R -- retrieve TCGA-BRCA expression and biospecimen data.
##
## Two queries: STAR counts with their colData, and the biospecimen aliquot
## table, used in 02 to check plate identifiers. Both are cached as .rds.
##
## GDCprepare() builds all six STAR assays at once (~3.6 GB), so files are
## prepared in chunks and only the `unstranded` counts are kept.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

suppressPackageStartupMessages(library(TCGAbiolinks))

## Files per GDCprepare() call; 150 keeps peak memory near 1 GB.
PREPARE_CHUNK_SIZE <- 150L

## GDCprepare() only adds an annotation column (such as
## paper_BRCA_Subtype_PAM50) to a chunk if some sample in that chunk has
## the annotation, so chunks can end up with different colData columns.
## cbind() needs identical columns. This adds any column a chunk lacks,
## filled with NA, so no annotation is lost when the chunks are combined.
alignColData <- function(pieces) {
    allCols <- unique(unlist(lapply(
        pieces, function(x) colnames(colData(x))
    )))

    ## Typed NAs, so filling does not coerce column types.
    proto <- list()

    for (cn in allCols) {
        for (p in pieces) {
            if (cn %in% colnames(colData(p))) {
                proto[[cn]] <- colData(p)[[cn]][0]
                break
            }
        }
    }

    lapply(pieces, function(p) {
        cd <- colData(p)

        for (cn in setdiff(allCols, colnames(cd))) {
            cd[[cn]] <- rep(proto[[cn]][NA_integer_], nrow(cd))
        }

        colData(p) <- cd[, allCols, drop = FALSE]
        p
    })
}

## ---- Expression ---------------------------------------------------------

if (file.exists(FILES$seRaw)) {
    say("01: ", basename(FILES$seRaw), " exists, skipping expression query")
} else {
    say("01: querying transcriptome profiling for ", CONFIG$project)

    queryExpr <- GDCquery(
        project       = CONFIG$project,
        data.category = "Transcriptome Profiling",
        data.type     = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
    )

    GDCdownload(queryExpr, directory = GDC_DIR)

    results <- getResults(queryExpr)
    chunks <- split(
        seq_len(nrow(results)),
        ceiling(seq_len(nrow(results)) / PREPARE_CHUNK_SIZE)
    )

    say(
        "01: preparing ", nrow(results), " files in ", length(chunks),
        " chunks of up to ", PREPARE_CHUNK_SIZE
    )

    pieces <- vector("list", length(chunks))

    for (i in seq_along(chunks)) {
        queryChunk <- queryExpr
        queryChunk$results[[1]] <- results[chunks[[i]], , drop = FALSE]

        ## Suppress GDCprepare()'s progress bar in the log.
        invisible(capture.output(
            seChunk <- suppressMessages(
                GDCprepare(queryChunk, directory = GDC_DIR)
            )
        ))

        ## Keep only the unstranded counts.
        if (!"unstranded" %in% assayNames(seChunk)) {
            stop(
                "Expected an 'unstranded' assay. Available: ",
                paste(assayNames(seChunk), collapse = ", ")
            )
        }

        assays(seChunk) <- assays(seChunk)["unstranded"]

        pieces[[i]] <- seChunk

        say(
            "01: chunk ", i, "/", length(chunks), " -> ",
            ncol(seChunk), " samples"
        )

        rm(seChunk, queryChunk)
        invisible(gc(verbose = FALSE))
    }

    nColsBefore <- vapply(
        pieces, function(x) ncol(colData(x)), integer(1)
    )

    pieces <- alignColData(pieces)

    if (length(unique(nColsBefore)) > 1L) {
        say(
            "01: colData columns per chunk ranged ", min(nColsBefore), "-",
            max(nColsBefore), "; aligned to ",
            ncol(colData(pieces[[1]])), " by union"
        )
    }

    ## All files share one gene annotation, so rowRanges match across chunks.
    seRaw <- do.call(cbind, pieces)

    rm(pieces)
    invisible(gc(verbose = FALSE))

    say(
        "01: prepared ", nrow(seRaw), " features x ",
        ncol(seRaw), " samples"
    )

    stopifnot(
        !anyDuplicated(colnames(seRaw)),
        identical(assayNames(seRaw), "unstranded")
    )

    saveCheckpoint(seRaw, FILES$seRaw)
    rm(seRaw)
    invisible(gc(verbose = FALSE))
}

## ---- Biospecimen (plate identifiers) ------------------------------------

if (file.exists(FILES$aliquot)) {
    say("01: ", basename(FILES$aliquot), " exists, skipping biospecimen query")
} else {
    say("01: querying biospecimen supplement (BCR Biotab)")

    queryBio <- GDCquery(
        project       = CONFIG$project,
        data.category = "Biospecimen",
        data.type     = "Biospecimen Supplement",
        data.format   = "BCR Biotab"
    )

    GDCdownload(queryBio, directory = GDC_DIR)
    bio <- suppressMessages(GDCprepare(queryBio, directory = GDC_DIR))

    ## Stop if the GDC renames the aliquot table.
    aliquotName <- "biospecimen_aliquot_brca"

    if (!aliquotName %in% names(bio)) {
        stop(
            "Expected '", aliquotName, "' in the biospecimen supplement. ",
            "Available: ", paste(names(bio), collapse = ", ")
        )
    }

    aliquot <- bio[[aliquotName]]

    say("01: aliquot table has ", nrow(aliquot), " rows")

    saveCheckpoint(aliquot, FILES$aliquot)
    rm(bio, aliquot)
    invisible(gc(verbose = FALSE))
}

say("01: done")
