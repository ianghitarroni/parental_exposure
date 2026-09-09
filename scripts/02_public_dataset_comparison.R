#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggwordcloud)
})

parse_args <- function(args) {
  required <- c(
    "--de-results",
    "--gse75613-table",
    "--gse169632-table",
    "--output"
  )
  missing_flags <- required[!required %in% args]
  if (length(missing_flags) > 0) {
    stop(
      paste0(
        "Missing required arguments: ", paste(missing_flags, collapse = ", "),
        "\nUsage: Rscript scripts/02_public_dataset_comparison.R ",
        "--de-results results/deseq2_significant_transcripts.csv ",
        "--gse75613-table data/table_sharma2016.csv ",
        "--gse169632-table data/translatome_embryo.csv ",
        "--output results/public_datasets"
      ),
      call. = FALSE
    )
  }

  values <- list()
  for (flag in required) {
    position <- match(flag, args)
    if (position == length(args) || startsWith(args[position + 1], "--")) {
      stop(paste("No value supplied for", flag), call. = FALSE)
    }
    values[[sub("^--", "", flag)]] <- args[position + 1]
  }
  values
}

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

assert_columns <- function(data, required, label) {
  missing_columns <- setdiff(required, colnames(data))
  if (length(missing_columns) > 0) {
    stop(
      label, " is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
assert_file(args$`de-results`, "Differential-results file")
assert_file(args$`gse75613-table`, "GSE75613 processed table")
assert_file(args$`gse169632-table`, "GSE169632 processed table")
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

de_results <- read.csv(
  args$`de-results`,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert_columns(de_results, c("symbol", "significant"), "Differential-results file")

de_symbols <- de_results %>%
  filter(significant, !is.na(symbol), symbol != "") %>%
  distinct(symbol) %>%
  pull(symbol)

if (length(de_symbols) == 0) {
  stop("No significant gene symbols were found.", call. = FALSE)
}

# GSE75613 comparison. The manuscript figure uses the seven columns below and
# excludes Blood, EpiMCA, ESC, Liver, Muscle and Testis from the heatmap.
gse75613 <- read.csv2(
  args$`gse75613-table`,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
gse75613_columns <- c(
  "Name", "Class", "EpiCA", "EpiCP", "SpCA", "SpCP",
  "Scyte", "Sptid1", "Sptid2"
)
assert_columns(gse75613, gse75613_columns, "GSE75613 processed table")

gse75613_overlap <- gse75613 %>%
  filter(Name %in% de_symbols, Class == "genes") %>%
  dplyr::select(all_of(gse75613_columns)) %>%
  distinct(Name, .keep_all = TRUE)

write.csv(
  gse75613_overlap,
  file.path(args$output, "gse75613_differential_transcript_overlap.csv"),
  row.names = FALSE
)

if (nrow(gse75613_overlap) >= 2) {
  gse75613_matrix <- gse75613_overlap %>%
    dplyr::select(-Name, -Class) %>%
    as.matrix()
  storage.mode(gse75613_matrix) <- "numeric"
  rownames(gse75613_matrix) <- gse75613_overlap$Name
  gse75613_matrix <- log2(gse75613_matrix + 1)

  png(
    file.path(args$output, "gse75613_expression_heatmap.png"),
    width = 1900,
    height = 3000,
    res = 250
  )
  pheatmap(
    gse75613_matrix,
    cutree_cols = 4,
    cutree_rows = 5,
    fontsize_row = 7.5,
    main = "Expression across reproductive cells and epididymal compartments"
  )
  dev.off()
}

# GSE169632 comparison with one-cell embryo total and ribosome-associated RNA.
gse169632 <- read.csv2(
  args$`gse169632-table`,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert_columns(
  gse169632,
  c("totalRNA", "translRNA", "TE", "symbol"),
  "GSE169632 processed table"
)

gse169632_overlap <- gse169632 %>%
  filter(is.finite(TE), TE > 0, symbol %in% de_symbols) %>%
  distinct(symbol, .keep_all = TRUE)

write.csv(
  gse169632_overlap,
  file.path(args$output, "gse169632_translatome_overlap.csv"),
  row.names = FALSE
)

cluster_autophagy <- c("Atg9a", "Ulk2", "Cib1", "Stk35", "Tbc1d5", "Epn1")
cluster_motility <- c("Tpgs2", "Lrrc8b")

if (nrow(gse169632_overlap) > 0) {
  wordcloud_data <- gse169632_overlap %>%
    mutate(
      functional_cluster = case_when(
        symbol %in% cluster_autophagy ~ "Autophagy-related cluster",
        symbol %in% cluster_motility ~ "Motility and osmo-adaptation cluster",
        TRUE ~ "Not assigned"
      )
    )

  wordcloud_plot <- ggplot(
    wordcloud_data,
    aes(label = symbol, size = TE, color = functional_cluster)
  ) +
    geom_text_wordcloud(
      rm_outside = TRUE,
      max_steps = 2,
      grid_size = 5,
      eccentricity = 0.9
    ) +
    scale_color_manual(
      values = c(
        "Not assigned" = "black",
        "Autophagy-related cluster" = "#D52954",
        "Motility and osmo-adaptation cluster" = "#515C9A"
      )
    ) +
    theme_void() +
    theme(legend.position = "none")

  ggsave(
    file.path(args$output, "gse169632_translatome_wordcloud.png"),
    wordcloud_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
}

public_dataset_check <- tibble(
  metric = c("gse75613_overlap", "gse169632_translatome_overlap"),
  expected = c(75L, 38L),
  observed = c(nrow(gse75613_overlap), nrow(gse169632_overlap)),
  status = if_else(expected == observed, "MATCH", "REVIEW")
)
write.csv(
  public_dataset_check,
  file.path(args$output, "public_dataset_result_check.csv"),
  row.names = FALSE
)

if (any(public_dataset_check$status == "REVIEW")) {
  warning(
    "Public-dataset overlaps do not fully match the manuscript summary. ",
    "Review public_dataset_result_check.csv and the processed input versions."
  )
}

capture.output(
  sessionInfo(),
  file = file.path(args$output, "sessionInfo_public_datasets.txt")
)
message("Public-dataset comparison completed: ", normalizePath(args$output))

