## run_all.R -- reproduce the TCGA-BRCA analysis in the BatchVaria
## application note.
##
##   Rscript run_all.R
##
## Each step can also be run on its own after sourcing R/00_setup.R, once the
## steps before it have run. Step 01 skips the download if data/raw/ exists.
## A full run from an existing download takes about 35 minutes. Step 10 runs
## last: it collects the numbers quoted in the manuscript from results/*.csv.

PROJECT_ROOT <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = TRUE),
    error = function(e) normalizePath(".", mustWork = TRUE)
)

STEPS <- c(
    "01_download_tcga.R",
    "02_sample_metadata.R",
    "03_preprocess.R",
    "04_build_object.R",
    "05_corrections.R",
    "06_profile_variance.R",
    "07_results_tables.R",
    "08_figures.R",
    "09_weighting_comparison.R",
    "10_manuscript_numbers.R"
)

for (step in STEPS) {
    message("\n=== ", step, " ===")
    source(file.path(PROJECT_ROOT, "R", step))
}

message("\n=== complete ===")
message("Tables in results/, figures in figures/")

## Record the package versions that produced these results.
writeLines(
    capture.output(sessionInfo()),
    file.path(PROJECT_ROOT, "results", "sessionInfo.txt")
)
