## 07_results_tables.R -- the tables the manuscript reports.
##
## Fractions, assay totals and the absolute contributions they imply
## (fraction x total, from full-precision fractions).

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

bv <- requireCheckpoint(FILES$bvProfiled, "R/06_profile_variance.R")

assaysPresent <- intersect(ASSAYS_ALL, assayNames(bv))

say("07: summarising ", paste(assaysPresent, collapse = ", "))

## ---- Total variance per assay ------------------------------------------

## The denominator for every fraction below.
av <- assayVariance(bv, assays = assaysPresent)

print(av)
writeTable(av, "07_assay_variance")

## ---- Fraction, total and absolute ----------------------------------------

## Valid as absolute contributions because fractions are pooled
## (CONFIG$anovaWeighting).
say("07: fraction, total and absolute variance (anova)")

vr <- varianceResults(bv, method = "anova")
vr <- vr[vr$assay %in% assaysPresent, c("assay", "term", "variance_fraction")]

headline <- merge(vr, av[, c("assay", "total_variance")], by = "assay")
headline$percent  <- 100 * headline$variance_fraction
headline$absolute <- headline$variance_fraction * headline$total_variance

## Assays in configured order, uncorrected first.
headline <- headline[
    order(headline$term, factor(headline$assay, levels = assaysPresent)),
    c("term", "assay", "percent", "total_variance", "absolute")
]

print(headline, digits = 6)
writeTable(headline, "07_anova_headline")

## ---- Where the biological fraction's rise comes from --------------------

## Split each change in the subtype fraction into numerator and
## denominator: f_new/f_old = (A_new/A_old) / (T_new/T_old). The identity is
## exact; `identity_check` is written out so it can be confirmed.
decomposeFractionChange <- function(headline, method) {
    rows <- headline[headline$term == CONFIG$bioVar, ]
    ref  <- rows[rows$assay == CONFIG$rawAssay, ]

    if (nrow(ref) != 1L) {
        say("07: no unique '", CONFIG$rawAssay, "' row to decompose against")
        return(NULL)
    }

    rows <- rows[rows$assay != CONFIG$rawAssay, ]

    out <- data.frame(
        method            = method,
        term              = CONFIG$bioVar,
        assay             = rows$assay,
        fraction_ratio    = rows$percent        / ref$percent,
        absolute_ratio    = rows$absolute       / ref$absolute,
        total_ratio       = rows$total_variance / ref$total_variance,
        stringsAsFactors  = FALSE
    )

    ## Percent changes, as quoted in the manuscript.
    out$fraction_change_pct <- 100 * (out$fraction_ratio - 1)
    out$absolute_change_pct <- 100 * (out$absolute_ratio - 1)
    out$total_change_pct    <- 100 * (out$total_ratio    - 1)

    out$identity_check <- out$fraction_ratio -
        out$absolute_ratio / out$total_ratio

    out
}

decomposition <- decomposeFractionChange(headline, "anova")

if (!is.null(decomposition)) {
    rownames(decomposition) <- NULL

    ## Fails if fractions and totals came from different objects.
    stopifnot(max(abs(decomposition$identity_check)) < 1e-8)

    say("07: decomposition of the change in the ", CONFIG$bioVar, " fraction")
    print(
        decomposition[, c(
            "method", "assay", "fraction_change_pct",
            "absolute_change_pct", "total_change_pct"
        )],
        digits = 4
    )

    writeTable(decomposition, "07_fraction_change_decomposition")
}

## ---- Basis retention ----------------------------------------------------

## Share of each assay's variance within the uncorrected assay's
## principal-axis subspace.
br <- basisRetention(bv, assays = assaysPresent)

say("07: basis retention")
print(br)
writeTable(br, "07_basis_retention")

say("07: done")
