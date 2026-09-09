#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

parse_args <- function(args) {
  required <- c("--counts", "--metadata", "--gtf", "--output")
  missing_flags <- required[!required %in% args]
  if (length(missing_flags) > 0) {
    stop(
      paste0(
        "Missing required arguments: ", paste(missing_flags, collapse = ", "),
        "\nUsage: Rscript scripts/01_deseq2_analysis.R ",
        "--counts 04_Counts/counts_matrix.txt ",
        "--metadata config/samples.csv ",
        "--gtf Ref_Genome_mm39/gencode.vM33.annotation.gtf ",
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
    dplyr::select(Geneid, all_of(sample_columns)) %>%
    column_to_rownames("Geneid") %>%
    as.matrix()

  storage.mode(matrix_data) <- "integer"
  matrix_data
}

extract_gtf_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;[[:space:]]*)", key, " \\\"([^\\\"]+)\\\"")
  stringr::str_match(attributes, pattern)[, 2]
}

read_gtf_annotations <- function(path) {
  if (!file.exists(path)) {
    stop("GTF annotation not found: ", path, call. = FALSE)
  }

  gtf <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE,
    col.names = c(
      "seqname", "source", "feature", "start", "end",
      "score", "strand", "frame", "attributes"
    )
  )

  gene_rows <- gtf %>% filter(feature == "gene")
  gene_type <- extract_gtf_attribute(gene_rows$attributes, "gene_type")
  missing_gene_type <- is.na(gene_type)
  gene_type[missing_gene_type] <- extract_gtf_attribute(
    gene_rows$attributes[missing_gene_type],
    "gene_biotype"
  )

  tibble(
    gene_id = extract_gtf_attribute(gene_rows$attributes, "gene_id"),
    gene_name = extract_gtf_attribute(gene_rows$attributes, "gene_name"),
    gene_type = gene_type
  ) %>%
    filter(!is.na(gene_id)) %>%
    mutate(gene_id_clean = str_remove(gene_id, "\\.[0-9]+$")) %>%
    distinct(gene_id_clean, .keep_all = TRUE) %>%
    dplyr::select(gene_id_clean, gene_name, gene_type)
}

normalize_sample_names <- function(names_vector) {
  names_vector %>%
    basename() %>%
    str_remove("\\.bam$") %>%
    str_remove("_Aligned.*$")
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
  mutate(condition = factor(condition, levels = c("Vehicle", "Cocaine"))) %>%
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

input_transcripts <- nrow(dds)
minimum_total_count <- 30L
keep <- rowSums(counts(dds)) >= minimum_total_count
dds <- dds[keep, ]
if (nrow(dds) == 0) {
  stop("No genes remain after the minimum total-count filter.", call. = FALSE)
}

dds <- DESeq(dds)
rld <- rlog(dds, blind = TRUE)

res <- results(
  dds,
  contrast = c("condition", "Cocaine", "Vehicle"),
  alpha = 0.1
)

gene_annotations <- read_gtf_annotations(args$gtf)

results_table <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_id_clean = str_remove(gene_id, "\\.[0-9]+$")) %>%
  left_join(gene_annotations, by = "gene_id_clean")

orgdb_symbols <- AnnotationDbi::mapIds(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = unique(results_table$gene_id_clean),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

results_table <- results_table %>%
  mutate(
    gencode_symbol = gene_name,
    symbol = unname(orgdb_symbols[gene_id_clean]),
    statistically_significant = !is.na(padj) &
      padj < 0.1 &
      abs(log2FoldChange) > 1,
    # The reference manuscript script maps Ensembl IDs with org.Mm.eg.db and
    # then applies na.omit(). Keep that reported set explicit while retaining
    # every statistically significant locus in the audit outputs.
    significant = statistically_significant & !is.na(symbol) & symbol != "",
    direction = case_when(
      significant & log2FoldChange > 1 ~ "Higher_in_cocaine",
      significant & log2FoldChange < -1 ~ "Lower_in_cocaine",
      TRUE ~ "Not_significant"
    ),
    biotype_class = case_when(
      gene_type == "protein_coding" ~ "Coding",
      !is.na(gene_type) ~ "Non-coding",
      TRUE ~ "Unknown"
    )
  ) %>%
  dplyr::select(-gene_name) %>%
  arrange(padj)

significant_loci_table <- results_table %>% filter(statistically_significant)
significant_table <- results_table %>% filter(significant)
excluded_unmapped_table <- results_table %>%
  filter(statistically_significant, !significant)

write.csv(
  results_table,
  file.path(args$output, "deseq2_all_transcripts.csv"),
  row.names = FALSE
)
write.csv(
  significant_table,
  file.path(args$output, "deseq2_significant_transcripts.csv"),
  row.names = FALSE
)
write.csv(
  significant_loci_table,
  file.path(args$output, "deseq2_significant_loci_before_symbol_mapping.csv"),
  row.names = FALSE
)
write.csv(
  excluded_unmapped_table,
  file.path(args$output, "excluded_unmapped_significant_loci.csv"),
  row.names = FALSE
)
writeLines(
  sort(unique(na.omit(significant_table$symbol))),
  file.path(args$output, "functional_gene_symbols.txt")
)
write.csv(
  as.data.frame(counts(dds, normalized = TRUE)) %>% rownames_to_column("gene_id"),
  file.path(args$output, "normalized_counts.csv"),
  row.names = FALSE
)

analysis_parameters <- tibble(
  parameter = c(
    "contrast", "minimum_total_count", "adjusted_p_value_threshold",
    "absolute_log2_fold_change_threshold", "multiple_testing",
    "exploratory_transformation", "rlog_blind"
  ),
  value = c(
    "Cocaine_vs_Vehicle", as.character(minimum_total_count), "0.1", "1",
    "Benjamini-Hochberg", "rlog", "TRUE"
  )
)
write.csv(
  analysis_parameters,
  file.path(args$output, "analysis_parameters.csv"),
  row.names = FALSE
)

direction_biotype_summary <- significant_table %>%
  count(biotype_class, direction, name = "transcripts") %>%
  arrange(biotype_class, direction)
write.csv(
  direction_biotype_summary,
  file.path(args$output, "direction_biotype_summary.csv"),
  row.names = FALSE
)

pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"), 1)
pca_data <- pca_data %>% rownames_to_column("sample")
write.csv(pca_data, file.path(args$output, "pca_coordinates.csv"), row.names = FALSE)

pca_plot <- ggplot(
  pca_data,
  aes(PC1, PC2, label = sample, color = condition, shape = condition)
) +
  geom_point(size = 4) +
  geom_text_repel(show.legend = FALSE, color = "black") +
  scale_color_manual(values = c("Vehicle" = "#2C7BB6", "Cocaine" = "#D7191C")) +
  labs(
    title = "PCA of mature spermatozoon RNA cargo",
    subtitle = "Regularized-logarithm transformed counts",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    color = "Condition",
    shape = "Condition"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(args$output, "pca_rlog.png"),
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

volcano_data <- results_table %>%
  mutate(
    plot_padj = case_when(
      is.na(padj) ~ NA_real_,
      padj <= 0 ~ .Machine$double.xmin,
      TRUE ~ padj
    ),
    label = if_else(
      significant,
      symbol,
      NA_character_
    )
  )

volcano_plot <- ggplot(
  volcano_data,
  aes(log2FoldChange, -log10(plot_padj), color = statistically_significant)
) +
  geom_point(alpha = 0.70, size = 1.7, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_text_repel(
    aes(label = label),
    max.overlaps = 20,
    size = 3,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "#C44E52")) +
  labs(
    title = "Cocaine versus vehicle",
    subtitle = "Differential representation in mature spermatozoon RNA cargo",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Meets thresholds"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(args$output, "volcano_plot.png"),
  volcano_plot,
  width = 7,
  height = 5,
  dpi = 300
)

if (nrow(direction_biotype_summary) > 0) {
  direction_biotype_plot <- direction_biotype_summary %>%
    mutate(
      signed_transcripts = if_else(
        direction == "Lower_in_cocaine",
        -transcripts,
        transcripts
      )
    ) %>%
    ggplot(aes(biotype_class, signed_transcripts, fill = direction)) +
    geom_col(width = 0.65) +
    coord_flip() +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(
      values = c(
        "Lower_in_cocaine" = "#2B83BA",
        "Higher_in_cocaine" = "#D01C8B"
      )
    ) +
    labs(
      title = "Direction and biotype of differential transcripts",
      x = NULL,
      y = "Number of transcripts",
      fill = NULL
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(args$output, "direction_by_biotype.png"),
    direction_biotype_plot,
    width = 7,
    height = 4.5,
    dpi = 300
  )
}

significant_gene_ids <- significant_table$gene_id
if (length(significant_gene_ids) >= 2) {
  heatmap_matrix <- assay(rld)[significant_gene_ids, , drop = FALSE]
  display_labels <- significant_table %>%
    arrange(match(gene_id, significant_gene_ids)) %>%
    transmute(label = coalesce(symbol, gene_id_clean)) %>%
    pull(label)
  rownames(heatmap_matrix) <- make.unique(display_labels)

  annotation <- data.frame(Condition = metadata$condition)
  rownames(annotation) <- rownames(metadata)

  png(
    file.path(args$output, "heatmap_differential_transcripts.png"),
    width = 1800,
    height = 3000,
    res = 250
  )
  pheatmap(
    heatmap_matrix,
    annotation_col = annotation,
    scale = "row",
    show_rownames = TRUE,
    fontsize_row = 7.5,
    cutree_rows = 2,
    main = "Differentially represented transcripts"
  )
  dev.off()
}

analysis_summary <- tibble(
  input_transcripts = input_transcripts,
  retained_after_count_filter = nrow(dds),
  statistically_significant_loci = sum(results_table$statistically_significant),
  excluded_unmapped_significant_loci = nrow(excluded_unmapped_table),
  significant_transcripts = sum(results_table$significant),
  higher_in_cocaine = sum(results_table$direction == "Higher_in_cocaine"),
  lower_in_cocaine = sum(results_table$direction == "Lower_in_cocaine"),
  coding_transcripts = sum(significant_table$biotype_class == "Coding"),
  noncoding_transcripts = sum(significant_table$biotype_class == "Non-coding"),
  unknown_biotype = sum(significant_table$biotype_class == "Unknown")
)
write.csv(
  analysis_summary,
  file.path(args$output, "analysis_summary.csv"),
  row.names = FALSE
)

manuscript_result_check <- tibble(
  metric = c(
    "significant_transcripts", "higher_in_cocaine", "lower_in_cocaine",
    "coding_transcripts", "noncoding_transcripts"
  ),
  expected = c(75L, 2L, 73L, 71L, 4L),
  observed = c(
    analysis_summary$significant_transcripts,
    analysis_summary$higher_in_cocaine,
    analysis_summary$lower_in_cocaine,
    analysis_summary$coding_transcripts,
    analysis_summary$noncoding_transcripts
  ),
  status = if_else(expected == observed, "MATCH", "REVIEW")
)
write.csv(
  manuscript_result_check,
  file.path(args$output, "manuscript_result_check.csv"),
  row.names = FALSE
)

if (any(manuscript_result_check$status == "REVIEW")) {
  warning(
    "Observed results do not fully match the manuscript summary. ",
    "Review results/manuscript_result_check.csv before updating figures or text."
  )
}

capture.output(sessionInfo(), file = file.path(args$output, "sessionInfo.txt"))
message("Analysis completed. Results written to: ", normalizePath(args$output))
