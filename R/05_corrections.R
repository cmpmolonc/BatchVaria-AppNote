## 05_corrections.R -- apply batch corrections.
##
## Four corrections: sva::ComBat() and limma::removeBatchEffect(), each
## naive (plate only) and subtype-protected (subtype named in `preserve`).

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

bv <- requireCheckpoint(FILES$bvRaw, "R/04_build_object.R")

set.seed(CONFIG$seed)

say("05: ", length(CORRECTIONS), " corrections on assay '", CONFIG$rawAssay, "'")

outcomes <- list()

for (nm in names(CORRECTIONS)) {
    spec <- CORRECTIONS[[nm]]

    say(
        "05: ", nm, " (method=", spec$method, ", preserve=",
        if (is.null(spec$preserve)) "none" else spec$preserve, ")"
    )

    noOpWarning <- FALSE

    ## Record a failed correction and continue with the others.
    result <- tryCatch(
        withCallingHandlers(
            runCorrection(
                bv,
                method       = spec$method,
                batch        = CONFIG$batchVar,
                preserve     = spec$preserve,
                assayName    = CONFIG$rawAssay,
                newAssayName = nm
            ),
            warning = function(w) {
                if (grepl("indistinguishable", conditionMessage(w))) {
                    noOpWarning <<- TRUE
                }

                say("05:   warning: ", conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) {
            say("05:   FAILED: ", conditionMessage(e))
            structure(conditionMessage(e), class = "correctionFailure")
        }
    )

    if (inherits(result, "correctionFailure")) {
        outcomes[[nm]] <- data.frame(
            assay = nm, method = spec$method,
            preserve = if (is.null(spec$preserve)) NA_character_ else spec$preserve,
            status = "failed", no_op = NA,
            message = as.character(result)
        )

        next
    }

    bv <- result

    outcomes[[nm]] <- data.frame(
        assay = nm, method = spec$method,
        preserve = if (is.null(spec$preserve)) NA_character_ else spec$preserve,
        status = "ok", no_op = noOpWarning,
        message = NA_character_
    )
}

say("05: assays now ", paste(assayNames(bv), collapse = ", "))

## ---- Correction ledger --------------------------------------------------

## Total variance before and after each correction, from the ledger.
ledger <- do.call(rbind, lapply(
    metadata(bv)$correction_history,
    function(e) {
        data.frame(
            assay_in           = e$assay_in,
            assay_out          = e$assay_out,
            method             = e$method,
            batch              = e$batch,
            preserve           = if (is.null(e$preserve)) NA_character_ else paste(e$preserve, collapse = "+"),
            no_op              = e$no_op,
            total_variance_in  = e$total_variance_in,
            total_variance_out = e$total_variance_out
        )
    }
))

print(ledger)

writeTable(do.call(rbind, outcomes), "05_correction_outcomes")
writeTable(ledger, "05_correction_ledger")

saveCheckpoint(bv, FILES$bvCorr)
say("05: done")
