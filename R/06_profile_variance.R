## 06_profile_variance.R -- decompose variance in every assay.
##
## anova engine: Type II sums of squares, pooled summarisation.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

bv <- requireCheckpoint(FILES$bvCorr, "R/05_corrections.R")

set.seed(CONFIG$seed)

## Only assays that 05 produced.
assaysToProfile <- intersect(ASSAYS_ALL, assayNames(bv))

say("06: profiling ", paste(assaysToProfile, collapse = ", "))

## ---- anova --------------------------------------------------------------

say(
    "06: anova, formula ", deparse1(CONFIG$anovaFormula),
    ", weighting '", CONFIG$anovaWeighting, "'"
)

t0 <- Sys.time()

## `weighting` is passed through to the engine (see CONFIG$anovaWeighting).
bv <- profileVariance(
    bv,
    formula   = CONFIG$anovaFormula,
    assays    = assaysToProfile,
    methods   = "anova",
    weighting = CONFIG$anovaWeighting
)

say("06: anova done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

## ---- Ledger check -------------------------------------------------------

## Warn if any assay is missing from the variance ledger.
history <- varianceHistory(bv)

recorded <- do.call(rbind, lapply(
    history,
    function(e) data.frame(assay = e$assay, method = e$method)
))

print(table(recorded$method, recorded$assay))

expected <- length(assaysToProfile)

if (nrow(recorded) < expected) {
    warning(
        "Expected ", expected, " ledger entries, found ", nrow(recorded),
        ". Some assays failed; see warnings above.",
        call. = FALSE
    )
}

writeTable(recorded, "06_ledger_entries")

saveCheckpoint(bv, FILES$bvProfiled)
say("06: done")
