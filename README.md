# Parental cocaine exposure — sperm RNA cargo analysis

Bioinformatic workflow developed for the undergraduate thesis **“Caracterización bioinformática del cargo de ARN del espermatozoide maduro asociado a la exposición paterna a cocaína”**.

**Author:** Ian Franco Ghitarroni  
**Institution:** Universidad Argentina de la Empresa (UADE)  
**Thesis director:** Dra. Betina González  
**Co-director:** Dr. Juan Antonio Bizzotto

## Scope

This repository contains the code used to analyze the **long-RNA/lncRNA-seq cargo of mature motile spermatozoa** recovered from the cauda epididymis and selected by swim-up. It compares two pooled biological replicates from vehicle-treated animals (`A2`, `A3`) with two pooled biological replicates from cocaine-treated animals (`B3`, `B5`).

The thesis scope is restricted to RNA-seq analysis. RRBS and multi-omic integration are not part of this repository.

## Workflow

1. Quality control with FastQC.
2. Adapter trimming and read filtering with fastp.
3. Alignment to the mouse reference genome GRCm39 using STAR.
4. Gene-level quantification with featureCounts and GENCODE vM33.
5. Low-count filtering and variance-stabilizing transformation.
6. Differential representation analysis with DESeq2.
7. PCA, volcano plot, heatmap and export of analysis-ready tables.
8. Downstream biological interpretation using Enrichr, STRING, GSE75613 and GSE169362.

The primary statistical contrast is **cocaine versus vehicle**. Following the reference manuscript, transcripts are classified as differentially represented when `padj < 0.1` and `|log2FoldChange| > 1`.

## Repository structure

```text
.
├── config/
│   └── samples.example.csv
├── scripts/
│   └── 01_deseq2_analysis.R
├── results/                 # generated locally; not versioned
├── .gitignore
└── README.md
```

## Input requirements

The R script expects:

- a featureCounts output file containing gene-level integer counts;
- a sample metadata CSV based on `config/samples.example.csv`;
- R 4.5.x or a compatible release;
- the packages listed in the installation command below.

Raw FASTQ, BAM, genome indexes and unpublished primary data are intentionally excluded because of file size, privacy, provenance and publication constraints.

## Installation

```r
install.packages(c("tidyverse", "pheatmap", "ggrepel"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "org.Mm.eg.db", "AnnotationDbi"))
```

## Execution

```bash
Rscript scripts/01_deseq2_analysis.R \
  --counts data/counts_matrix.txt \
  --metadata config/samples.csv \
  --output results
```

The output directory includes the complete DESeq2 table, the filtered significant-transcript table, normalized counts, PCA coordinates, a PCA figure, a volcano plot, a top-transcript heatmap and `sessionInfo.txt`.

## Interpretation constraints

- The mature spermatozoon is not treated as a conventionally transcriptionally active cell. Results are therefore described as changes in **RNA cargo representation**, not as active gene expression.
- `n = 2` per condition corresponds to two independent pooled biological replicates; each pool was generated from six animals. It must not be interpreted as twelve independent RNA-seq observations.
- The workflow identifies statistical associations in the sperm RNA cargo. It does not, by itself, establish causality, paternal transmission or a specific offspring phenotype.
- Detection of a transcript in one-cell embryo ribosome-associated fractions is compatible with early translational availability, but does not prove exclusive paternal origin.
- Biological interpretation should remain aligned with the reference manuscript and distinguish spermatogenic from epididymal contributions.

## Data availability

This repository does not redistribute raw sequencing data or third-party datasets. Public reference datasets must be obtained from their original repositories under their respective terms. Analysis-ready inputs may be added only after confirming authorization, publication status and applicable repository limits.

## Reproducibility and provenance

The pipeline exports R session information with every run. For a frozen computational environment, create an `renv` lockfile from the validated analysis workstation before final publication.

## Citation

Until the associated thesis and manuscript receive their final bibliographic records, cite this repository using its GitHub URL, author, title and accessed commit hash.

## License and reuse

No open-source license has been assigned yet. Consequently, the code is publicly visible but no reuse rights are granted beyond those provided by applicable law. A license should be selected only after confirming the requirements of UADE, the research group and the associated manuscript.