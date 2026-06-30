#!/usr/bin/env Rscript

# Annotate a consensus peak BED against a reference GTF with ChIPseeker.
#
# For each peak we report its genomic feature class (promoter / UTR / exon /
# intron / intergenic, etc.), the nearest gene, and the distance to that gene's
# TSS. A TxDb is built on the fly from the GTF.
#
# Two outputs are written:
#   --bed : a BED6+ file (0-based, IGV/bedtools friendly) with the nearest
#           gene as the name and the feature class + distance-to-TSS appended.
#   --tsv : an optional, header-bearing table with the same annotations plus
#           gene-level coordinates that don't fit cleanly in a BED.

suppressPackageStartupMessages({
  library(optparse)
  library(ChIPseeker)
  library(txdbmaker)
})

# columns emitted to the .tsv (annotatePeak produces more; we keep a stable set)
ANNOTATED_TSV_COLUMNS <- c(
  "seqnames", "start", "end", "strand",
  "annotation", "geneId", "distanceToTSS",
  "geneStart", "geneEnd", "geneStrand"
)

parse_arguments <- function() {
  option_list <- list(
    make_option("--peaks", type = "character",
                help = "Input consensus peak .bed file for one mark."),
    make_option("--gtf", type = "character",
                help = "Reference annotation .gtf file."),
    make_option("--bed", type = "character",
                help = "Path to write the annotated .bed file (BED6+)."),
    make_option("--tsv", type = "character", default = NULL,
                help = "Optional path to write an annotated .tsv with a header."),
    make_option("--upstream", type = "integer", default = 3000,
                help = "Bases upstream of the TSS for the promoter region [default %default]."),
    make_option("--downstream", type = "integer", default = 3000,
                help = "Bases downstream of the TSS for the promoter region [default %default].")
  )
  arguments <- parse_args(OptionParser(option_list = option_list))

  required <- c("peaks", "gtf", "bed")
  missing <- required[vapply(required, function(name) is.null(arguments[[name]]), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste0("--", missing, collapse = ", "))
  }
  arguments
}

# Convert annotatePeak's data frame (1-based, GRanges-style) into a BED6+ frame:
# chrom, start (0-based), end, name (nearest gene), score, strand, then the
# feature class and distance-to-TSS as extra columns.
build_bed <- function(annotated_df) {
  strand <- as.character(annotated_df$strand)
  gene_id <- as.character(annotated_df$geneId)
  data.frame(
    chrom         = as.character(annotated_df$seqnames),
    start         = annotated_df$start - 1L,
    end           = annotated_df$end,
    name          = ifelse(is.na(gene_id), ".", gene_id),
    score         = ".",
    strand        = ifelse(strand %in% c("+", "-"), strand, "."),
    annotation    = as.character(annotated_df$annotation),
    distanceToTSS = annotated_df$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

write_bed <- function(bed, output_path) {
  write.table(bed, file = output_path, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_tsv <- function(table, output_path) {
  write.table(table, file = output_path, sep = "\t",
              quote = FALSE, row.names = FALSE)
}

ensure_dir <- function(path) {
  output_dir <- dirname(path)
  if (nzchar(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

main <- function() {
  arguments <- parse_arguments()

  ensure_dir(arguments$bed)
  if (!is.null(arguments$tsv)) {
    ensure_dir(arguments$tsv)
  }

  # an empty consensus (no peaks survived the threshold for this mark) still
  # needs outputs so downstream rules have a consistent, present file set
  if (length(readLines(arguments$peaks, warn = FALSE)) == 0) {
    file.create(arguments$bed)  # an empty .bed is valid
    if (!is.null(arguments$tsv)) {
      empty <- data.frame(matrix(nrow = 0, ncol = length(ANNOTATED_TSV_COLUMNS)))
      colnames(empty) <- ANNOTATED_TSV_COLUMNS
      write_tsv(empty, arguments$tsv)
    }
    message("No consensus peaks; wrote empty annotation for ", arguments$bed)
    return(invisible())
  }

  txdb <- txdbmaker::makeTxDbFromGFF(arguments$gtf, format = "gtf")
  peaks <- readPeakFile(arguments$peaks)
  annotated <- annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-arguments$upstream, arguments$downstream),
    level = "gene",
    verbose = FALSE
  )
  annotated_df <- as.data.frame(annotated)

  write_bed(build_bed(annotated_df), arguments$bed)
  if (!is.null(arguments$tsv)) {
    keep <- intersect(ANNOTATED_TSV_COLUMNS, colnames(annotated_df))
    write_tsv(annotated_df[, keep, drop = FALSE], arguments$tsv)
  }
  message("Annotated ", nrow(annotated_df), " peaks -> ", arguments$bed)
}

main()
