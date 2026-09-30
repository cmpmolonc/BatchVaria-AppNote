## 10_manuscript_numbers.R -- the numbers quoted in the manuscript.
##
## Reads results/*.csv and writes results/10_manuscript_numbers.md, each
## block naming the file it came from. Runs last; needs only the tables from
## 02-09, not the profiled object.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

## ---- Reading -------------------------------------------------------------

## Stop if a table is missing.
readResult <- function(name) {
    path <- file.path(PATHS$results, paste0(name, ".csv"))

    if (!file.exists(path)) {
        stop(
            "results/", name, ".csv is missing. Run the step that writes it ",
            "before 10.", call. = FALSE
        )
    }

    read.csv(path, stringsAsFactors = FALSE)
}

design     <- readResult("02_design_summary")
dataSource <- readResult("02_data_source")
sampleFlow <- readResult("02_sample_flow")
filtering  <- readResult("03_feature_filtering")
outcomes   <- readResult("05_correction_outcomes")
totals     <- readResult("07_assay_variance")
retention  <- readResult("07_basis_retention")
anovaHead  <- readResult("07_anova_headline")
decomp     <- readResult("07_fraction_change_decomposition")
weighting  <- readResult("09_weighting_absolute")
nullLevel  <- readResult("09_plate_null_level")

## ---- Assay order and labels ---------------------------------------------

assayOrder <- totals$assay[order(match(totals$assay, ASSAYS_ALL))]

LABELS <- c(
    raw              = "uncorrected",
    combat_naive     = "ComBat naive",
    combat_protected = "ComBat protected",
    limma_naive      = "limma naive",
    limma_protected  = "limma protected"
)

lab <- function(a) ifelse(a %in% names(LABELS), LABELS[a], a)

## ---- Formatting ---------------------------------------------------------

## Four significant figures; values that are zero up to floating point are
## printed in scientific notation.
fmt <- function(x, digits = 4) {
    small <- is.finite(x) & x != 0 & abs(x) < 1e-4

    ifelse(
        small,
        formatC(x, format = "e", digits = 1),
        formatC(x, format = "fg", digits = digits, big.mark = "", flag = "#")
    )
}

mdTable <- function(df, digits = 4) {
    body <- df

    ## Integers print as integers.
    for (j in seq_along(body)) {
        if (!is.numeric(body[[j]])) next

        v <- body[[j]]
        isCount <- all(is.na(v) | (is.finite(v) & v == round(v)))
        body[[j]] <- if (isCount) formatC(v, format = "d") else fmt(v, digits)
    }

    c(
        paste0("| ", paste(colnames(body), collapse = " | "), " |"),
        paste0("|", paste(rep("---", ncol(body)), collapse = "|"), "|"),
        apply(body, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
    )
}

pick <- function(df, term, assay, column) {
    df[[column]][df$term == term & df$assay == assay]
}

out <- character(0)
add <- function(...) out <<- c(out, ...)

## ---- Header -------------------------------------------------------------

add(
    "# Numbers for the manuscript",
    "",
    paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"),
           " by `R/10_manuscript_numbers.R` from the tables in `results/`."),
    "",
    "Every number below is read from a table in `results/`.",
    "",
    paste0("Variance engine: `anova`, Type II sums of squares, ",
           "`weighting = \"", CONFIG$anovaWeighting, "\"`."),
    ""
)

## ---- Dataset ------------------------------------------------------------

add(
    "## Dataset",
    "",
    "Source: `02_design_summary.csv`, `03_feature_filtering.csv`.",
    "",
    paste0("- **", design$n_samples, " samples**, ", design$n_plates,
           " sequencing plates, ", design$n_subtypes, " PAM50 subtypes."),
    paste0("- Smallest retained plate: ", design$min_plate_size, " samples."),
    paste0("- Features: ", filtering$n_features[filtering$step == "input"],
           " measured, **",
           filtering$n_features[filtering$step == "variance_filter"],
           " retained** after `filterByExpr`."),
    paste0("- Plate/subtype association: Cramer's V = **", design$cramers_v,
           "**, chi-squared p = ", format(design$chisq_p, digits = 3), "."),
    ""
)

## ---- Dataset provenance -------------------------------------------------

## Data source and sample flow.

release <- dataSource$value[dataSource$field == "gdc_data_release"]

add(
    "## Dataset provenance",
    "",
    "Source: `02_data_source.csv`, `02_sample_flow.csv`.",
    "",
    mdTable(dataSource),
    "",
    "Cohort reduction, in the order applied:",
    "",
    mdTable(sampleFlow),
    ""
)

if (is.na(release) || !nzchar(release)) {
    add(
        paste0("The GDC data release is not recorded; data retrieved on ",
               dataSource$value[dataSource$field == "retrieved"], "."),
        ""
    )
}

## ---- Corrections --------------------------------------------------------

retention <- retention[order(match(retention$assay, ASSAYS_ALL)), ]

add(
    "## Corrections applied",
    "",
    "Source: `05_correction_outcomes.csv`, `07_basis_retention.csv`.",
    "",
    ## Counts read from the outcomes table.
    paste0(
        nrow(outcomes), " corrections attempted; ",
        sum(outcomes$status == "ok"), " completed, ",
        sum(outcomes$no_op), " were no-ops.",
        if (any(outcomes$status != "ok")) paste0(
            " **Failed: ",
            paste(outcomes$assay[outcomes$status != "ok"], collapse = ", "),
            ".**"
        ) else "",
        if (any(outcomes$no_op)) paste0(
            " **No-op: ",
            paste(outcomes$assay[outcomes$no_op], collapse = ", "),
            ".**"
        ) else ""
    ),
    "",
    mdTable(data.frame(
        assay     = lab(retention$assay),
        retention = retention$retention,
        off_basis = retention$off_basis
    ), digits = 6),
    ""
)

## ---- The headline table -------------------------------------------------

anovaHead <- anovaHead[order(anovaHead$term, match(anovaHead$assay, ASSAYS_ALL)), ]

bioRows <- anovaHead[anovaHead$term == CONFIG$bioVar, ]

add(
    "## Subtype variance",
    "",
    "Source: `07_anova_headline.csv`.",
    "",
    paste0("Subtype variance under Type II sums of squares, as a fraction ",
           "and as an absolute quantity."),
    "",
    mdTable(data.frame(
        assay            = lab(bioRows$assay),
        `fraction (%)`   = bioRows$percent,
        `assay total`    = bioRows$total_variance,
        `absolute`       = bioRows$absolute,
        check.names      = FALSE
    ), digits = 7),
    ""
)

## Assays whose absolute subtype contribution equals the uncorrected value.
bioAbs <- setNames(bioRows$absolute, bioRows$assay)
refAbs <- bioAbs[[CONFIG$rawAssay]]
sameAbs <- names(bioAbs)[abs(bioAbs - refAbs) / refAbs < 1e-8]

if (length(sameAbs) > 1L) {
    sameRows <- bioRows[bioRows$assay %in% sameAbs, ]

    add(
        paste0("Assays with the same absolute subtype variance as the ",
               "uncorrected assay (", fmt(refAbs, 7), "):"),
        "",
        mdTable(data.frame(
            assay          = lab(sameRows$assay),
            `fraction (%)` = sameRows$percent,
            `assay total`  = sameRows$total_variance,
            absolute       = sameRows$absolute,
            check.names    = FALSE
        ), digits = 7),
        ""
    )
}

## ---- Where the fraction's rise comes from -------------------------------

decomp <- decomp[order(match(decomp$assay, ASSAYS_ALL)), ]

add(
    "## Change in the subtype fraction",
    "",
    "Source: `07_fraction_change_decomposition.csv`.",
    "",
    "Percent change from uncorrected; `f_new/f_old = (A_new/A_old) / (T_new/T_old)`.",
    "",
    mdTable(data.frame(
        assay             = lab(decomp$assay),
        `fraction (%)`    = decomp$fraction_change_pct,
        `absolute (%)`    = decomp$absolute_change_pct,
        `total (%)`       = decomp$total_change_pct,
        check.names       = FALSE
    ), digits = 4),
    "",
    paste0("The identity closes to ", format(max(abs(decomp$identity_check)),
                                             digits = 2),
           "."),
    ""
)

## ---- shared -------------------------------------------------------------

sharedRows <- anovaHead[anovaHead$term == "shared", ]
plateRows  <- anovaHead[anovaHead$term == CONFIG$batchVar, ]

add(
    "## Signed shared remainder",
    "",
    "Source: `07_anova_headline.csv`.",
    "",
    ## Plate column included because the limma identity below uses it.
    mdTable(data.frame(
        assay          = lab(sharedRows$assay),
        `shared (%)`   = sharedRows$percent,
        `shared (abs)` = sharedRows$absolute,
        `plate (abs)`  = plateRows$absolute[match(sharedRows$assay, plateRows$assay)],
        check.names    = FALSE
    ), digits = 5),
    ""
)

## Under naive removeBatchEffect, `shared` = -plate exactly; checked here.
if ("limma_naive" %in% sharedRows$assay) {
    lnShared <- pick(sharedRows, "shared", "limma_naive", "absolute")
    lnPlate  <- pick(plateRows, CONFIG$batchVar, "limma_naive", "absolute")

    relDiff <- abs(lnShared + lnPlate) / abs(lnPlate)

    identityHolds <- relDiff < 1e-8

    ## Report how closely they cancel.
    sigFigs <- if (relDiff > 0) floor(-log10(relDiff)) else Inf

    add(
        paste0("Under limma naive, `shared` = ", fmt(lnShared, 12),
               " and `plate` = ", fmt(lnPlate, 12),
               ". They sum to ", formatC(lnShared + lnPlate, format = "e",
                                         digits = 1),
               ", i.e. they are equal and opposite to ",
               if (is.finite(sigFigs)) paste(sigFigs, "significant figures")
               else "the last bit",
               ". This holds by construction for naive `removeBatchEffect`;",
               " the ComBat value is empirical."),
        ""
    )

    stopifnot(identityHolds)
}

## ---- Is the corrected plate term actually small? ------------------------

nullLevel <- nullLevel[order(match(nullLevel$assay, ASSAYS_ALL)), ]

add(
    "## Plate term against its noise-only level",
    "",
    "Source: `09_plate_null_level.csv`.",
    "",
    paste0("Null level `df_plate/df_resid` = ",
           fmt(nullLevel$null_ratio[1], 6), " (",
           nullLevel$df_plate[1], "/", nullLevel$df_resid[1], ")."),
    "",
    mdTable(data.frame(
        assay              = lab(nullLevel$assay),
        `plate / residual` = nullLevel$plate_over_residual,
        `x null level`     = nullLevel$x_null,
        check.names        = FALSE
    ), digits = 4),
    ""
)

## ---- What averaging fractions costs -------------------------------------

wRaw <- weighting[weighting$assay == CONFIG$rawAssay, ]

bioErr   <- wRaw$error_pct[wRaw$term == CONFIG$bioVar]
plateErr <- wRaw$error_pct[wRaw$term == CONFIG$batchVar]

weighting <- weighting[order(weighting$term, match(weighting$assay, ASSAYS_ALL)), ]
wBio <- weighting[weighting$term == CONFIG$bioVar, ]

add(
    "## Feature-averaged and pooled summaries",
    "",
    "Source: `09_weighting_absolute.csv`.",
    "",
    paste0("On uncorrected data the feature summary **understates subtype by ",
           fmt(abs(bioErr), 3), "%** and **overstates plate by ",
           fmt(abs(plateErr), 3), "%**."),
    "",
    "Subtype variance under each summary:",
    "",
    mdTable(data.frame(
        assay       = lab(wBio$assay),
        feature     = wBio$feature,
        pooled      = wBio$pooled,
        `error (%)` = wBio$error_pct,
        check.names = FALSE
    ), digits = 6),
    ""
)

## ---- Figures ------------------------------------------------------------

add(
    "## Figures",
    "",
    "| figure | drawn from | script |",
    "|---|---|---|",
    "| `figure1_main` A | schematic of the workflow; no data | `R/08_figures.R` |",
    "| `figure1_main` B | `07_anova_headline.csv`: fraction, total, absolute and `shared` | `R/08_figures.R` |",
    "| `figure1_main` C | `07_fraction_change_decomposition.csv` | `R/08_figures.R` |",
    "",
    "Panels B and C are Type II sums of squares under pooled weighting.",
    ""
)

## ---- Claims this analysis does not support ------------------------------

add(
    "## Not supported by these results",
    "",
    "- That naive and protected limma differ in absolute subtype variance.",
    "  They are identical, and provably so.",
    "- That either naive correction leaves a plate effect above chance. Both",
    "  are below the null level for a term of this size.",
    "- That ComBat and limma independently confirm `shared` as a diagnostic.",
    "  One is a measurement; the other is an identity.",
    ""
)

path <- file.path(PATHS$results, "10_manuscript_numbers.md")
writeLines(out, path)

say("10: wrote ", path)
say("10: done")
