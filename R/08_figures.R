## 08_figures.R -- Figure 1.
##
##   A  BatchVaria workflow (schematic, no data)
##   B  subtype fraction, total variance, absolute contribution and `shared`
##   C  percent change in fraction, absolute contribution and total variance

source(file.path(PROJECT_ROOT, "R", "00_setup.R"))

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

bv <- requireCheckpoint(FILES$bvProfiled, "R/06_profile_variance.R")

assaysPresent <- intersect(ASSAYS_ALL, assayNames(bv))

## ---- Presentation constants --------------------------------------------

## Colours: blue subtype, orange plate, aqua `shared`, grey residual.
TERM_COLOURS <- c(
    plate_id                 = "#eb6834",
    paper_BRCA_Subtype_PAM50 = "#2a78d6",
    shared                   = "#1baf7a",
    residual                 = "#b9b7ae"
)

METHOD_COLOURS <- c(
    Uncorrected = "#2a78d6",
    ComBat      = "#eb6834",
    limma       = "#1baf7a"
)

ASSAY_LABELS <- c(
    raw              = "Uncorrected",
    combat_naive     = "ComBat\nnaive",
    combat_protected = "ComBat\nprotected",
    limma_naive      = "limma\nnaive",
    limma_protected  = "limma\nprotected"
)

ASSAY_METHOD <- c(
    raw              = "Uncorrected",
    combat_naive     = "ComBat",
    combat_protected = "ComBat",
    limma_naive      = "limma",
    limma_protected  = "limma"
)

theme_bv <- function(base_size = 8) {
    theme_minimal(base_size = base_size) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
            axis.text = element_text(colour = "grey25"),
            axis.title = element_text(colour = "grey20"),
            strip.text = element_text(face = "bold", colour = "grey15"),
            legend.title = element_blank(),
            legend.key.size = unit(0.8, "lines"),
            plot.title = element_text(face = "bold", size = rel(1.0)),
            plot.subtitle = element_text(colour = "grey35", size = rel(0.85)),
            plot.caption = element_text(colour = "grey40", size = 7, hjust = 0)
        )
}

decorate <- function(df, assayCol = "assay") {
    df$assay_label <- factor(
        ASSAY_LABELS[as.character(df[[assayCol]])],
        levels = ASSAY_LABELS[assaysPresent]
    )

    df$method_group <- factor(
        ASSAY_METHOD[as.character(df[[assayCol]])],
        levels = names(METHOD_COLOURS)
    )

    df
}

savePlot <- function(plot, name, width, height) {
    ## cairo_pdf, so UTF-8 characters render in the PDF.
    devices <- list(pdf = grDevices::cairo_pdf, png = "png")

    for (ext in names(devices)) {
        withCallingHandlers(
            ggsave(
                file.path(PATHS$figures, paste0(name, ".", ext)),
                plot,
                width = width, height = height, dpi = 300,
                device = devices[[ext]]
            ),
            ## Stop rather than silently substitute a glyph.
            warning = function(w) {
                if (grepl("mbcsToSbcs|no font could be found|invalid.*glyph",
                          conditionMessage(w))) {
                    stop(
                        "glyph substitution while writing ", name, ".", ext,
                        ": ", conditionMessage(w), call. = FALSE
                    )
                }
            }
        )
    }

    say("08: wrote figures/", name, ".{pdf,png}")
}

## ---- Type sizes -----------------------------------------------------------

## Drawn at printed size (18 cm, double column). geom_text() sizes are in mm;
## PT converts from points. Nothing is below 6.5 pt.
PT <- 1 / ggplot2::.pt
TXT_BODY  <- 7.0 * PT
TXT_SMALL <- 6.5 * PT
TXT_TITLE <- 8.5

FIG_WIDTH_IN  <- 18 / 2.54
FIG_HEIGHT_IN <- 16 / 2.54

## ---- Data ---------------------------------------------------------------

av <- assayVariance(bv, assays = assaysPresent)

## Absolute contribution = pooled fraction x assay total, at full precision.
buildFrame <- function(method) {
    vr <- varianceResults(bv, method = method)
    vr <- vr[vr$assay %in% assaysPresent, c("assay", "term", "variance_fraction")]

    df <- merge(vr, av[, c("assay", "total_variance")], by = "assay")

    df$percent <- 100 * df$variance_fraction
    df$absolute <- df$variance_fraction * df$total_variance
    df$assay <- factor(df$assay, levels = assaysPresent)

    decorate(df)
}

anovaDf <- buildFrame("anova")

## =========================================================================
## Figure 1
## =========================================================================

## ---- Panel A: the workflow ----------------------------------------------

## Schematic stacks on an absolute-variance axis: after correction the plate
## block shrinks and the subtype block keeps its height.

## Block heights. `bio` is shared by both stacks so the subtype blocks match.
bio <- 1.45
plateRaw <- 1.25
plateCorrected <- 0.10
remainderRaw <- 3.90
remainderCorrected <- 3.40

barBase <- 1.30

leftX <- c(1.45, 3.35)
rightX <- c(12.75, 14.65)

## Cumulative tops, so a change to any block moves everything above it.
stackTops <- function(base, heights) base + cumsum(heights)

rawTops <- stackTops(barBase, c(remainderRaw, bio, plateRaw))
corTops <- stackTops(barBase, c(remainderCorrected, bio, plateCorrected))

schematic <- data.frame(
    xmin = c(rep(leftX[1], 3), rep(rightX[1], 3)),
    xmax = c(rep(leftX[2], 3), rep(rightX[2], 3)),
    ymin = c(barBase, rawTops[1:2], barBase, corTops[1:2]),
    ymax = c(rawTops, corTops),
    term = rep(c("residual", CONFIG$bioVar, CONFIG$batchVar), 2),
    stringsAsFactors = FALSE
)

schematic$term <- factor(
    schematic$term,
    levels = c("residual", CONFIG$bioVar, CONFIG$batchVar)
)

## The two subtype blocks must be the same height.
bioHeights <- with(
    schematic[schematic$term == CONFIG$bioVar, ], ymax - ymin
)

stopifnot(diff(range(bioHeights)) < 1e-9)
stopifnot(max(corTops) < max(rawTops))

## Block labels; the thin corrected plate block gets a leader line instead.
blockLabels <- data.frame(
    x = c(mean(leftX), mean(leftX), mean(leftX), mean(rightX), mean(rightX)),
    y = c(
        barBase + remainderRaw / 2,
        rawTops[1] + bio / 2,
        rawTops[2] + plateRaw / 2,
        barBase + remainderCorrected / 2,
        corTops[1] + bio / 2
    ),
    label = c(
        "Residual", "Subtype", "Plate",
        "Residual", "Subtype"
    ),
    stringsAsFactors = FALSE
)

boxX <- c(5.15, 11.15)

## The box's top edge is level with the top of the uncorrected stack; its
## contents, arrows and loop move with it (BOX_DY). A positive curvature bends
## the upward loop away from the box.
BOX_TOP <- 9.35
BOX_DY <- rawTops[3] - BOX_TOP
ARROW_Y <- 4.35 + BOX_DY
LOOP_CURVATURE <- 0.55
boxMidX <- mean(boxX)

panelA <- ggplot() +
    coord_cartesian(xlim = c(0, 16), ylim = c(0, 10), expand = FALSE) +
    ## ---- the two stacks
    geom_rect(
        data = schematic,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = term),
        colour = "grey20", linewidth = 0.35
    ) +
    geom_text(
        data = blockLabels,
        aes(x, y, label = label),
        colour = "white", size = TXT_BODY
    ) +
    ## Leader line for the thin corrected plate block.
    annotate(
        "segment",
        x = rightX[2] + 0.12, xend = rightX[2] + 0.55,
        y = corTops[2] + plateCorrected / 2, yend = corTops[3] + 0.55,
        linewidth = 0.3, colour = "grey45"
    ) +
    annotate(
        "text", x = rightX[2] + 0.62, y = corTops[3] + 0.68,
        label = "Plate", hjust = 0, size = TXT_BODY, colour = "grey30"
    ) +
    ## ---- the absolute axis
    annotate(
        "segment", x = 0.72, xend = 0.72, y = barBase, yend = 9.15,
        arrow = arrow(length = unit(0.22, "cm")),
        linewidth = 0.5, colour = "grey30"
    ) +
    annotate(
        "text", x = 0.34, y = 5.1, angle = 90,
        label = "Absolute variance", size = TXT_BODY, colour = "grey25"
    ) +
    annotate(
        "text", x = mean(leftX), y = 0.72,
        label = "Uncorrected", size = TXT_BODY, colour = "grey20"
    ) +
    annotate(
        "text", x = mean(rightX), y = 0.72,
        label = "Corrected", size = TXT_BODY, colour = "grey20"
    ) +
    ## ---- the BatchVaria box
    annotate(
        "rect", xmin = boxX[1], xmax = boxX[2], ymin = (3.85 + BOX_DY), ymax = (9.35 + BOX_DY),
        fill = "grey96", colour = "grey30", linewidth = 0.55
    ) +
    annotate(
        "text", x = boxMidX, y = (8.97 + BOX_DY),
        label = "BatchVaria", fontface = "bold", size = 8.5 * PT, colour = "grey10"
    ) +
    annotate(
        "text", x = boxMidX, y = (8.45 + BOX_DY),
        label = "Provenance-tracked variance decomposition",
        size = TXT_SMALL, colour = "grey30"
    ) +
    annotate(
        "segment", x = boxX[1] + 0.2, xend = boxX[2] - 0.2,
        y = (8.05 + BOX_DY), yend = (8.05 + BOX_DY), linewidth = 0.35, colour = "grey55"
    ) +
    annotate(
        "text", x = boxMidX, y = (7.66 + BOX_DY),
        label = "profileVariance()", size = TXT_BODY, colour = "grey15"
    ) +
    annotate(
        "text", x = boxMidX, y = (7.13 + BOX_DY),
        label = "fraction, absolute variance, shared",
        size = TXT_SMALL, colour = "grey40", fontface = "italic"
    ) +
    annotate(
        "segment", x = boxX[1] + 0.2, xend = boxX[2] - 0.2,
        y = (6.73 + BOX_DY), yend = (6.73 + BOX_DY), linewidth = 0.35, colour = "grey55"
    ) +
    annotate(
        "text", x = boxMidX, y = (6.34 + BOX_DY),
        label = "runCorrection()", size = TXT_BODY, colour = "grey15"
    ) +
    annotate(
        "text", x = boxMidX, y = (5.81 + BOX_DY),
        label = "ComBat or limma; naive or subtype-protected",
        size = TXT_SMALL, colour = "grey40", fontface = "italic"
    ) +
    annotate(
        "segment", x = boxX[1] + 0.2, xend = boxX[2] - 0.2,
        y = (5.41 + BOX_DY), yend = (5.41 + BOX_DY), linewidth = 0.35, colour = "grey55"
    ) +
    annotate(
        "text", x = boxMidX, y = (5.02 + BOX_DY),
        label = "re-profile and evaluate", size = TXT_BODY, colour = "grey15"
    ) +
    annotate(
        "text", x = boxMidX, y = (4.49 + BOX_DY),
        label = "assayVariance(), basisRetention()",
        size = TXT_SMALL, colour = "grey40", fontface = "italic"
    ) +
    ## ---- arrows in and out (below the loop, so they do not cross it)
    annotate(
        "segment", x = leftX[2] + 0.15, xend = boxX[1] - 0.15,
        y = ARROW_Y, yend = ARROW_Y,
        arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
        linewidth = 0.6, colour = "grey35"
    ) +
    annotate(
        "segment", x = boxX[2] + 0.15, xend = rightX[1] - 0.15,
        y = ARROW_Y, yend = ARROW_Y,
        arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
        linewidth = 0.6, colour = "grey35"
    ) +
    ## ---- the evaluation loop: up the box's right-hand side, from
    ## re-profiling back to profiling
    annotate(
        "curve", x = boxX[2] + 0.08, xend = boxX[2] + 0.08,
        y = (5.02 + BOX_DY), yend = (7.66 + BOX_DY),
        curvature = LOOP_CURVATURE,
        arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
        linewidth = 0.55, colour = "grey35"
    ) +
    annotate(
        "text", x = boxMidX, y = (3.25 + BOX_DY),
        label = "Variance evaluation workflow",
        size = TXT_BODY, colour = "grey25"
    ) +
    scale_fill_manual(values = TERM_COLOURS, guide = "none") +
    theme_void() +
    theme(
        plot.title = element_text(
            colour = "grey20", size = TXT_TITLE, face = "bold", hjust = 0,
            margin = margin(b = 4)
        ),
        plot.margin = margin(4, 4, 2, 4)
    )

## ---- Panel B: the table -------------------------------------------------

## Drawn as a table within the figure; values from the anova results.
tableTerms <- c(CONFIG$bioVar, "shared")

tableWide <- reshape(
    anovaDf[anovaDf$term %in% tableTerms, c("assay", "term", "percent", "absolute")],
    idvar = "assay", timevar = "term", direction = "wide"
)

tableWide <- tableWide[match(assaysPresent, tableWide$assay), ]

tableWide$total <- av$total_variance[match(tableWide$assay, av$assay)]

bioPct <- tableWide[[paste0("percent.", CONFIG$bioVar)]]
bioAbs <- tableWide[[paste0("absolute.", CONFIG$bioVar)]]
## `shared` is shown as a percentage, as in the text.
sharedAbs <- tableWide[["absolute.shared"]]
sharedPct <- tableWide[["percent.shared"]]

stopifnot(
    !anyNA(bioPct), !anyNA(bioAbs), !anyNA(sharedAbs), !anyNA(sharedPct),
    !anyNA(tableWide$total),
    ## Same sign as a fraction and as an absolute.
    identical(sign(sharedPct), sign(sharedAbs))
)

## Rows whose absolute contribution equals the uncorrected value (bolded).
isInvariant <- abs(bioAbs - bioAbs[tableWide$assay == CONFIG$rawAssay]) /
    bioAbs[tableWide$assay == CONFIG$rawAssay] < 1e-8

nRow <- nrow(tableWide)

## Geometry. Row 0 is the header; rows run downward.
colX <- c(assay = 0.015, pct = 0.455, total = 0.665, abs = 0.845, shared = 0.995)
rowY <- seq(from = nRow, by = -1, length.out = nRow)
headY <- nRow + 0.72

cell <- function(x, y, label, hjust, fontface = "plain", colour = "grey20",
                 size = TXT_BODY, vjust = 0.5) {
    data.frame(
        x = x, y = y, label = label, hjust = hjust, vjust = vjust,
        fontface = fontface, colour = colour, size = size,
        stringsAsFactors = FALSE
    )
}

cells <- rbind(
    cell(colX[["assay"]],  headY, "Assay",                   0, "bold", "grey10", vjust = 0),
    cell(colX[["pct"]],    headY, "Subtype\nfraction (%)",   1, "bold", "grey10", vjust = 0),
    cell(colX[["total"]],  headY, "Total\nvariance",         1, "bold", "grey10", vjust = 0),
    cell(colX[["abs"]],    headY, "Subtype\nabsolute",       1, "bold", "grey10", vjust = 0),
    cell(colX[["shared"]], headY, "Shared\n(%)",             1, "bold", "grey10", vjust = 0),

    cell(colX[["assay"]], rowY, sub("\n", " ", ASSAY_LABELS[tableWide$assay]), 0),
    cell(colX[["pct"]],   rowY, sprintf("%.2f", bioPct), 1),
    ## One decimal, so 49 999.6 does not read as a rounded 50 000.
    cell(
        colX[["total"]], rowY,
        formatC(tableWide$total, format = "f", digits = 1, big.mark = " "), 1
    ),

    ## Bold blue where equal to the uncorrected value.
    cell(
        colX[["abs"]], rowY, sprintf("%.1f", bioAbs), 1,
        fontface = ifelse(isInvariant, "bold", "plain"),
        colour = ifelse(isInvariant, TERM_COLOURS[[CONFIG$bioVar]], "grey20")
    ),

    ## `shared`, signed.
    cell(colX[["shared"]], rowY, sprintf("%+.2f", sharedPct), 1)
)

panelB <- ggplot() +
    ## y padding keeps row spacing even when the panel is stretched.
    coord_cartesian(
        xlim = c(-0.02, 1.02), ylim = c(-4.4, headY + 2.2), expand = FALSE
    ) +
    geom_text(
        data = cells,
        aes(x, y, label = label, hjust = hjust, vjust = vjust,
            fontface = fontface, colour = colour, size = size),
        lineheight = 0.9
    ) +
    scale_colour_identity() +
    scale_size_identity() +
    annotate(
        "segment", x = -0.01, xend = 1.01, y = headY - 0.22, yend = headY - 0.22,
        linewidth = 0.5, colour = "grey35"
    ) +
    annotate(
        "segment", x = -0.01, xend = 1.01, y = 0.42, yend = 0.42,
        linewidth = 0.4, colour = "grey60"
    ) +
    theme_void() +
    theme(
        plot.title = element_text(
            colour = "grey20", size = TXT_TITLE, face = "bold", hjust = 0,
            margin = margin(b = 6)
        ),
        plot.margin = margin(4, 6, 4, 6)
    )

## ---- Panel C: where the rise comes from ---------------------------------

## From 07's fraction-change decomposition.
decompFile <- file.path(PATHS$results, "07_fraction_change_decomposition.csv")

if (!file.exists(decompFile)) {
    stop("results/07_fraction_change_decomposition.csv is missing; run 07.",
         call. = FALSE)
}

decompAll <- read.csv(decompFile, stringsAsFactors = FALSE)
decompDf <- decompAll[decompAll$method == "anova", ]

stopifnot(nrow(decompDf) > 0)

MEASURE_LEVELS <- c(
    "Fraction",
    "Absolute (numerator)",
    "Total (denominator)"
)

changeDf <- rbind(
    data.frame(assay = decompDf$assay, measure = MEASURE_LEVELS[1],
               value = decompDf$fraction_change_pct),
    data.frame(assay = decompDf$assay, measure = MEASURE_LEVELS[2],
               value = decompDf$absolute_change_pct),
    data.frame(assay = decompDf$assay, measure = MEASURE_LEVELS[3],
               value = decompDf$total_change_pct)
)

changeDf$measure <- factor(changeDf$measure, levels = MEASURE_LEVELS)
changeDf <- decorate(changeDf)

## Hollow bars for the fraction, filled for the absolute, grey for the total.
MEASURE_FILL <- setNames(
    c("white", TERM_COLOURS[[CONFIG$bioVar]], "#b9b7ae"), MEASURE_LEVELS
)

MEASURE_LINE <- setNames(
    c(TERM_COLOURS[[CONFIG$bioVar]], TERM_COLOURS[[CONFIG$bioVar]], "#8f8d86"),
    MEASURE_LEVELS
)

## Absolute changes to two decimals; exact zeros print as "0".
isAbsolute <- changeDf$measure == MEASURE_LEVELS[2]
isZero <- abs(changeDf$value) < 1e-9

changeLabel <- ifelse(
    isZero, "0",
    ifelse(isAbsolute, sprintf("%+.2f", changeDf$value),
           sprintf("%+.1f", changeDf$value))
)

## Tick marks for exact-zero bars (the middle bar, so no dodge offset).
stopifnot(all(changeDf$measure[isZero] == MEASURE_LEVELS[2]))

zeroTicks <- data.frame(
    x = as.numeric(droplevels(changeDf$assay_label)[isZero])
)

panelC <- ggplot(changeDf, aes(assay_label, value, fill = measure, colour = measure)) +
    geom_col(
        position = position_dodge(width = 0.78),
        width = 0.70, linewidth = 0.45
    ) +
    geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey30") +
    geom_text(
        aes(label = changeLabel, vjust = ifelse(value >= 0 | isZero, -0.55, 1.45)),
        position = position_dodge(width = 0.78),
        size = TXT_SMALL, colour = "grey25", show.legend = FALSE
    ) +
    geom_segment(
        data = zeroTicks,
        aes(x = x - 0.12, xend = x + 0.12, y = 0, yend = 0),
        inherit.aes = FALSE,
        colour = TERM_COLOURS[[CONFIG$bioVar]], linewidth = 1.4
    ) +
    scale_fill_manual(values = MEASURE_FILL) +
    scale_colour_manual(values = MEASURE_LINE) +
    scale_y_continuous(expand = expansion(mult = c(0.14, 0.16))) +
    labs(x = NULL, y = "% change vs uncorrected") +
    theme_bv() +
    theme(
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.margin = margin(b = -4),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        plot.title = element_text(
            colour = "grey20", size = TXT_TITLE, face = "bold", hjust = 0
        )
    )

## ---- Assemble -----------------------------------------------------------

## A across the top; B and C below.
figure1 <- (panelA / (panelB | panelC)) +
    plot_layout(heights = c(0.82, 1)) +
    plot_annotation(
        tag_levels = "A",
        theme = theme_bv() +
            theme(plot.margin = margin(6, 6, 6, 6))
    ) &
    theme(plot.tag = element_text(face = "bold", size = 10))

savePlot(figure1, "figure1_main", width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN)

say("08: done")
