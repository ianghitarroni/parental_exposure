#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
})

parse_args <- function(args) {
  required <- c("--counts", "--metadata", "--output")
  missing_flags <- required[!required %in% args]
  if (length(missing_flags) > 0) {
    stop(
      paste0(
        "Missing required arguments: ", paste(missing_flags, collapse = ", "),
        "\nUsage: Rscript scripts/01_deseq2_analysis.R ",
        "--counts data/counts_matrix.txt ",
        "--metadata config/samples.csv ",
        "--output results"
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

read_featurecounts <- function(path) {
  if (!file.exists(path)) {
    stop("Counts file not found: ", path, call. = FALSE)
  }

  counts_raw <- read.delim(
    path,
    header = TRUE,
    comment.char = "#",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  missing_annotation <- setdiff(annotation_columns, colnames(counts_raw))
  if (length(missing_annotation) > 0) {
    stop(
      "The input does not look like a featureCounts table. Missing columns: ",
      paste(missing_annotation, collapse = ", "),
      call. = FALSE
    )
  }

  sample_columns <- setdiff(colnames(counts_raw), annotation_columns)
  if (length(sample_columns) < 2) {
    stop("At least two sample count columns are required.", call. = FALSE)
  }

  matrix_data <- counts_raw %>%
    select(Geneid, all_of(sample_columns)) %>%
    column_to_rownames("Geneid") %>%
    as.matrix()

  storage.mode(matrix_data) <- "integer"
  matrix_data
}

normalize_sample_names <- function(names_vector) {
  names_vector %>%
    basename() %>%
    str_remove("\\.bam$") %>%
    str_remove("\\.Aligned.*$")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(args$metadata, stringsAsFactors = FALSE, check.names = FALSE)
required_metadata <- c("sample", "condition")
if (!all(required_metadata %in% colnames(metadata))) {
  stop("Metadata must contain columns: sample, condition", call. = FALSE)
}
if (anyDuplicated(metadata$sample)) {
  stop("Metadata contains duplicated sample identifiers.", call. = FALSE)
}

count_matrix <- read_featurecounts(args$counts)
colnames(count_matrix) <- normalize_sample_names(colnames(count_matrix))

missing_in_counts <- setdiff(metadata$sample, colnames(count_matrix))
extra_in_counts <- setdiff(colnames(count_matrix), metadata$sample)
if (length(missing_in_counts) > 0 || length(extra_in_counts) > 0) {
  stop(
    paste0(
      "Counts and metadata do not match. ",
      "Missing in counts: ", paste(missing_in_counts, collapse = ", "), "; ",
      "Missing in metadata: ", paste(extra_in_counts, collapse = ", ")
    ),
    call. = FALSE
  )
}

metadata <- metadata %>%
  mutate(
    condition = factor(condition, levels = c("Vehicle", "Cocaine"))
  ) %>%
  arrange(match(sample, colnames(count_matrix))) %>%
  column_to_rownames("sample")

if (any(is.na(metadata$condition))) {
  stop("Condition values must be exactly Vehicle or Cocaine.", call. = FALSE)
}
if (any(table(metadata$condition) < 2)) {
  warning("At least one condition has fewer than two samples; inference will be fragile.")
}

count_matrix <- count_matrix[rownames(count_matrix) != "", , drop = FALSE]
count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]

if (any(count_matrix < 0, na.rm = TRUE) || anyNA(count_matrix)) {
  stop("Count matrix must contain non-negative integer values without NA.", call. = FALSE)
}

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
if (nrow(dds) == 0) {
  stop("No genes remain after the minimum total-count filter.", call. = FALSE)
}

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

res <- results(
  dds,
  contrast = c("condition", "Cocaine", "Vehicle"),
  alpha = 0.1
)

results_table <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_id_clean = str_remove(gene_id, "\\.[0-9]+$"),
    significant = !is.na(padj) & padj < 0.1 & abs(log2FoldChange) > 1,
    direction = case_when(
      significant & log2FoldChange > 1 ~ "Higher_in_cocaine",
      significant & log2FoldChange < -1 ~ "Lower_in_cocaine",
      TRUE ~ "Not_significant"
    )
  ) %>%
  arrange(padj)

write.csv(
  results_table,
  file.path(args$output, "deseq2_all_transcripts.csv"),
  row.names = FALSE
)
write.csv(
  filter(results_table, significant),
  file.path(args$output, "deseq2_significant_transcripts.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(counts(dds, normalized = TRUE)) %>% rownames_to_column("gene_id"),
  file.path(args$output, "normalized_counts.csv"),
  row.names = FALSE
)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"), 1)
pca_data <- pca_data %>% rownames_to_column("sample")
write.csv(pca_data, file.path(args$output, "pca_coordinates.csv"), row.names = FALSE)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, label = sample, shape = condition)) +
  geom_point(size = 4) +
  geom_text_repel(show.legend = FALSE) +
  labs(
    title = "PCA of mature sperm RNA cargo",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    shape = "Condition"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(args$output, "pca.png"),
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

volcano_data <- results_table %>%
  mutate(
    plot_pvalue = if_else(is.na(pvalue) | pvalue <= 0, .Machine$double.xmin, pvalue),
    label = if_else(significant, gene_id_clean, NA_character_)
  )

volcano_plot <- ggplot(
  volcano_data,
  aes(log2FoldChange, -log10(plot_pvalue), shape = significant)
) +
  geom_point(alpha = 0.65, size = 1.6, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_text_repel(
    aes(label = label),
    max.overlaps = 15,
    size = 3,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "Cocaine versus vehicle",
    subtitle = "Differential representation in mature sperm RNA cargo",
    x = "log2 fold change",
    y = "-log10 p-value",
    shape = "Meets thresholds"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(args$output, "volcano_plot.png"),
  volcano_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ranked_genes <- results_table %>%
  filter(!is.na(padj)) %>%
  slice_head(n = min(50, n())) %>%
  pull(gene_id)

if (length(ranked_genes) >= 2) {
  heatmap_matrix <- assay(vsd)[ranked_genes, , drop = FALSE]
  annotation <- data.frame(Condition = metadata$condition)
  rownames(annotation) <- rownames(metadata)

  png(
    file.path(args$output, "heatmap_top_transcripts.png"),
    width = 1800,
    height = 2200,
    res = 250
  )
  pheatmap(
    heatmap_matrix,
    annotation_col = annotation,
    scale = "row",
    show_rownames = FALSE,
    main = "Top transcripts by adjusted p-value"
  )
  dev.off()
}

summary_table <- results_table %>%
  summarise(
    tested_transcripts = n(),
    significant_transcripts = sum(significant),
    higher_in_cocaine = sum(direction == "Higher_in_cocaine"),
    lower_in_cocaine = sum(direction == "Lower_in_cocaine")
  )
write.csv(summary_table, file.path(args$output, "analysis_summary.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(args$output, "sessionInfo.txt"))
message("Analysis completed. Results written to: ", normalizePath(args$output))
