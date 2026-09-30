## 02_sample_metadata.R -- attach sequencing plate, and define the sample set.
##
## The sequencing plate is read from the RNA-seq aliquot barcode
## (TCGA-A8-A09E-01A-11R-A00Z-07: plate A00Z, characters 22-25) and checked
## against the biospecimen aliquot table. All sample selection happens here,
## before filtering and normalisation.

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))


seRaw <- requireCheckpoint(FILES$seRaw, "R/01_download_tcga.R")
aliquot <- requireCheckpoint(FILES$aliquot, "R/01_download_tcga.R")

## Sample counts at each step are written to results/02_sample_flow.csv.
nStart <- ncol(seRaw)
say("02: starting from ", nrow(seRaw), " features x ", nStart, " samples")

## ---- Primary tumours only ----------------------------------------------

seRaw <- seRaw[, seRaw$shortLetterCode == CONFIG$sampleType]
nPrimary <- ncol(seRaw)
say("02: ", nPrimary, " samples after restricting to ", CONFIG$sampleType)

## ---- Sequencing plate -------------------------------------------------

## Each expression column is one RNA aliquot; its barcode names the plate it
## was sequenced on. (A sample has several aliquots, often on different
## plates, so the 16-character sample barcode does not identify the plate.)
barcode <- colnames(seRaw)
colData(seRaw)[[CONFIG$batchVar]] <- substr(barcode, 22, 25)
colData(seRaw)$patient <- substr(barcode, 1, 12)

## Every aliquot found in the biospecimen table must carry the same plate.
tablePlate <- aliquot$plate_id[match(barcode, aliquot$bcr_aliquot_barcode)]
found <- !is.na(tablePlate)

say(
    "02: plate read from barcode; ", sum(found), " of ", length(barcode),
    " aliquots found in the biospecimen table"
)

stopifnot(identical(
    tablePlate[found],
    seRaw[[CONFIG$batchVar]][found]
))

## ---- One aliquot per patient ---------------------------------------------

## Record replicate aliquots before keeping one per patient.
patient <- seRaw$patient
dupPatients <- unique(patient[duplicated(patient)])

replicates <- data.frame(
    patient = patient[patient %in% dupPatients],
    barcode = barcode[patient %in% dupPatients],
    plate   = seRaw[[CONFIG$batchVar]][patient %in% dupPatients],
    stringsAsFactors = FALSE
)
replicates <- replicates[order(replicates$patient, replicates$barcode), ]

ord <- order(barcode, decreasing = identical(CONFIG$replicateRule, "latest"))
keep <- barcode[ord][!duplicated(patient[ord])]
replicates$kept <- replicates$barcode %in% keep

seRaw <- seRaw[, barcode %in% keep]
nReplicate <- nPrimary - ncol(seRaw)

say(
    "02: ", length(dupPatients), " patients with more than one aliquot; ",
    "removed ", nReplicate, " by rule '", CONFIG$replicateRule, "'"
)

stopifnot(!anyDuplicated(seRaw$patient))

writeTable(replicates, "02_replicate_aliquots")

## ---- Complete annotation ------------------------------------------------

hasBio <- !is.na(seRaw[[CONFIG$bioVar]])
hasBatch <- !is.na(seRaw[[CONFIG$batchVar]])

nNoBio <- sum(!hasBio)
nNoBatch <- sum(hasBio & !hasBatch)

say(
    "02: dropping ", nNoBio, " samples without ", CONFIG$bioVar,
    " and ", nNoBatch, " further without ", CONFIG$batchVar
)

seRaw <- seRaw[, hasBio & hasBatch]

## ---- Minimum plate size -------------------------------------------------

## Plates too small for ComBat to estimate a batch variance.
plateCounts <- table(seRaw[[CONFIG$batchVar]])
smallPlates <- names(plateCounts)[plateCounts < CONFIG$minPlateSize]

nSmallPlates <- length(smallPlates)
nSmallPlateSamples <- sum(plateCounts[smallPlates])

say(
    "02: ", nSmallPlates, " of ", length(plateCounts),
    " plates have < ", CONFIG$minPlateSize, " samples (",
    nSmallPlateSamples, " samples); removing"
)

seRaw <- seRaw[, !seRaw[[CONFIG$batchVar]] %in% smallPlates]

## ---- Factors ------------------------------------------------------------

## Model both variables as factors throughout.
colData(seRaw)[[CONFIG$batchVar]] <-
    factor(as.character(seRaw[[CONFIG$batchVar]]))
colData(seRaw)[[CONFIG$bioVar]] <-
    factor(as.character(seRaw[[CONFIG$bioVar]]))

## ---- Design diagnostics -------------------------------------------------

## Plate-subtype association (Cramer's V), reported in Section 3.
designTab <- table(
    plate = seRaw[[CONFIG$batchVar]],
    subtype = seRaw[[CONFIG$bioVar]]
)

chi <- suppressWarnings(chisq.test(designTab))
cramersV <- sqrt(
    as.numeric(chi$statistic) /
        (sum(designTab) * (min(dim(designTab)) - 1))
)

say(
    "02: final set ", nrow(seRaw), " features x ", ncol(seRaw), " samples, ",
    nlevels(seRaw[[CONFIG$batchVar]]), " plates, ",
    nlevels(seRaw[[CONFIG$bioVar]]), " subtypes"
)
say("02: plate-subtype association Cramer's V = ", round(cramersV, 4))

writeTable(
    data.frame(
        n_samples   = ncol(seRaw),
        n_features  = nrow(seRaw),
        n_plates    = nlevels(seRaw[[CONFIG$batchVar]]),
        n_subtypes  = nlevels(seRaw[[CONFIG$bioVar]]),
        min_plate_size = CONFIG$minPlateSize,
        cramers_v   = round(cramersV, 4),
        chisq_p     = signif(chi$p.value, 3)
    ),
    "02_design_summary"
)

writeTable(
    as.data.frame(designTab, responseName = "n"),
    "02_plate_by_subtype"
)

## ---- Provenance ---------------------------------------------------------

## Data source. The GDC data release cannot be recovered by the pipeline and
## is left NA; the retrieval date is recorded.
writeTable(
    data.frame(
        field = c(
            "project", "data_category", "data_type", "workflow_type",
            "sample_type_code", "batch_variable", "batch_source",
            "biological_variable", "biological_source",
            "retrieved", "gdc_data_release"
        ),
        value = c(
            CONFIG$project,
            "Transcriptome Profiling",
            "Gene Expression Quantification",
            "STAR - Counts",
            CONFIG$sampleType,
            CONFIG$batchVar,
            "RNA-seq aliquot barcode (characters 22-25), checked against GDC Biospecimen Supplement (BCR Biotab)",
            CONFIG$bioVar,
            "TCGAbiolinks colData; PAM50 calls from Berger et al. 2018 (Cancer Cell)",
            format(file.mtime(FILES$seRaw), "%Y-%m-%d"),
            NA_character_
        ),
        stringsAsFactors = FALSE
    ),
    "02_data_source"
)

## Samples remaining after each selection step, in order.
writeTable(
    data.frame(
        step = c(
            "downloaded",
            paste0("sample type == ", CONFIG$sampleType),
            paste0("one aliquot per patient (", CONFIG$replicateRule, ")"),
            paste0("has ", CONFIG$bioVar),
            paste0("has ", CONFIG$batchVar),
            paste0("plate size >= ", CONFIG$minPlateSize)
        ),
        n_samples = c(
            nStart,
            nPrimary,
            nPrimary - nReplicate,
            nPrimary - nReplicate - nNoBio,
            nPrimary - nReplicate - nNoBio - nNoBatch,
            ncol(seRaw)
        ),
        n_removed = c(
            NA_integer_,
            nStart - nPrimary,
            nReplicate,
            nNoBio,
            nNoBatch,
            nSmallPlateSamples
        ),
        stringsAsFactors = FALSE
    ),
    "02_sample_flow"
)

saveCheckpoint(seRaw, FILES$seAnnot)
say("02: done")
