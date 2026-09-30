## 09_weighting_comparison.R -- feature-averaged vs pooled summarisation,
## and the plate term against its noise-only level.
##
## Compares the two anova summaries ("feature" and "pooled") on every assay
## and term, as quoted in the Discussion, and asserts that pooled subtype and
## residual contributions are unchanged by either removeBatchEffect
## correction. Uses the corrected object from 05.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

suppressPackageStartupMessages({
    library(matrixStats)
})

bv <- requireCheckpoint(FILES$bvCorr, "R/05_corrections.R")

set.seed(CONFIG$seed)

assaysPresent <- intersect(ASSAYS_ALL, assayNames(bv))

## Call the engine directly so both summaries come from the same fit.
anovaEngine <- get(".varianceEngines", envir = asNamespace("BatchVaria"))[["anova"]]

sampleData <- as.data.frame(colData(bv))

say("09: comparing 'feature' and 'pooled' summaries on ",
    length(assaysPresent), " assays")

rows <- list()

for (assayName in assaysPresent) {
    m <- assay(bv, assayName)

    ## Same denominator as 07.
    totalVariance <- sum(rowVars(m))

    for (w in c("feature", "pooled")) {
        res <- anovaEngine(
            assayMatrix = m,
            formula     = CONFIG$anovaFormula,
            sampleData  = sampleData,
            weighting   = w
        )

        rows[[length(rows) + 1L]] <- data.frame(
            assay     = assayName,
            weighting = w,
            term      = res$term,
            fraction  = res$variance_fraction,
            total     = totalVariance,
            absolute  = res$variance_fraction * totalVariance,
            stringsAsFactors = FALSE
        )
    }

    say("09: ", assayName)
}

comparison <- do.call(rbind, rows)
writeTable(comparison, "09_weighting_comparison")

wide <- reshape(
    comparison[, c("assay", "term", "weighting", "absolute")],
    idvar = c("assay", "term"), timevar = "weighting", direction = "wide"
)

colnames(wide) <- sub("^absolute\\.", "", colnames(wide))

wide$error_pct <- 100 * (wide$feature - wide$pooled) / wide$pooled
wide$assay <- factor(wide$assay, levels = assaysPresent)
wide <- wide[order(wide$term, wide$assay), ]

say("09: absolute variance under each summary")
print(wide, row.names = FALSE, digits = 6)

writeTable(wide, "09_weighting_absolute")

## ---- The identity that settles which summary is right -------------------

## Pooled subtype and residual contributions must be identical in the
## uncorrected and both removeBatchEffect assays.
limmaAssays <- intersect(
    c("limma_naive", "limma_protected"), as.character(assaysPresent)
)

if (length(limmaAssays) && CONFIG$rawAssay %in% assaysPresent) {
    invariant <- c(CONFIG$bioVar, "residual")

    pooled <- comparison[comparison$weighting == "pooled", ]

    for (tm in intersect(invariant, unique(pooled$term))) {
        refValue <- pooled$absolute[
            pooled$term == tm & pooled$assay == CONFIG$rawAssay
        ]

        for (la in limmaAssays) {
            got <- pooled$absolute[pooled$term == tm & pooled$assay == la]

            say(
                "09: ", tm, " -- raw ", format(refValue, digits = 10),
                " vs ", la, " ", format(got, digits = 10)
            )

            ## Relative tolerance for floating-point error.
            stopifnot(abs(got - refValue) / refValue < 1e-8)
        }
    }

    say("09: removeBatchEffect invariance holds under 'pooled'")
}

## ---- Is the corrected batch term small, or merely small-looking? --------

## Under no plate effect, E[SS_plate] / E[RSS] = df_plate / df_resid; the
## plate-to-residual ratio is reported relative to that level.
designMatrix <- model.matrix(
    CONFIG$anovaFormula, data = sampleData
)

dfPlate <- nlevels(sampleData[[CONFIG$batchVar]]) - 1L
dfResid <- nrow(sampleData) - qr(designMatrix)$rank
nullRatio <- dfPlate / dfResid

say(
    "09: df_plate ", dfPlate, ", df_resid ", dfResid,
    ", null E[plate/residual] ", format(nullRatio, digits = 6)
)

pooled <- comparison[comparison$weighting == "pooled", ]

ratioOf <- function(assayName) {
    pooled$absolute[pooled$term == CONFIG$batchVar & pooled$assay == assayName] /
        pooled$absolute[pooled$term == "residual" & pooled$assay == assayName]
}

nullLevels <- do.call(rbind, lapply(assaysPresent, function(a) {
    r <- ratioOf(a)

    say(
        "09: ", a, " plate/residual ", format(r, digits = 6),
        " = ", format(r / nullRatio, digits = 4), "x the null level"
    )

    data.frame(
        assay              = a,
        plate_over_residual = r,
        null_ratio         = nullRatio,
        x_null             = r / nullRatio,
        df_plate           = dfPlate,
        df_resid           = dfResid,
        stringsAsFactors   = FALSE
    )
}))

## Written to results/09_plate_null_level.csv.
writeTable(nullLevels, "09_plate_null_level")

## The uncorrected data carry a plate effect above the null level.
stopifnot(ratioOf(CONFIG$rawAssay) > nullRatio)

## Every correction leaves it below.
correctedAssays <- setdiff(assaysPresent, CONFIG$rawAssay)
stopifnot(vapply(correctedAssays, function(a) ratioOf(a) < nullRatio, logical(1)))

say("09: plate above chance before correction, below chance after, in every assay")

say("09: done")
