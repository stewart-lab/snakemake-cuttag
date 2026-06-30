#!/usr/bin/env Rscript

# Annotate a consensus peak BED against a reference GTF with ChIPseeker.
#
# For each peak we report its genomic feature class (promoter / UTR / exon /
# intron / intergenic, etc.), the nearest gene, and the distance to that gene's
# TSS. A TxDb is built on the fly from the GTF.

suppressPackageStartupMessages({
  library(optparse)
  library(ChIPseeker)
  library(txdbmaker)
})

# the stable subset of columns we emit (annotatePeak produces more)
ANNOTATED_COLUMNS <- c(
  "seqnames", "start", "end", "strand",
  "annotation", "geneId", "distanceToTSS"
)

parse_arguments <- function() {
  option_list <- list(
    make_option("--peaks", type = "character",
                help = "Consensus peak .bed file for one mark."),
    make_option("--gtf", type = "character",
                help = "Reference annotation .gtf file."),
    make_option("--output", type = "character",
                help = "Path to write the annotated .tsv file."),
    make_option("--upstream", type = "integer", default = 3000,
                help = "Bases upstream of the TSS for the promoter region [default %default]."),
    make_option("--downstream", type = "integer", default = 3000,
                help = "Bases downstream of the TSS for the promoter region [default %default].")
  )
  arguments <- parse_args(OptionParser(option_list = option_list))

  required <- c("peaks", "gtf", "output")
  missing <- required[vapply(required, function(name) is.null(arguments[[name]]), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste0("--", missing, collapse = ", "))
  }
  arguments
}

write_tsv <- function(table, output_path) {
  write.table(table, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

main <- function() {
  arguments <- parse_arguments()

  output_dir <- dirname(arguments$output)
  if (nzchar(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # an empty consensus (no peaks survived the threshold for this mark) still
  # needs a (header-only) output so downstream rules have a consistent schema
  if (length(readLines(arguments$peaks, warn = FALSE)) == 0) {
    empty <- data.frame(matrix(nrow = 0, ncol = length(ANNOTATED_COLUMNS)))
    colnames(empty) <- ANNOTATED_COLUMNS
    write_tsv(empty, arguments$output)
    message("No consensus peaks; wrote header-only annotation: ", arguments$output)
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
  keep <- intersect(ANNOTATED_COLUMNS, colnames(annotated_df))
  write_tsv(annotated_df[, keep, drop = FALSE], arguments$output)
  message("Annotated ", nrow(annotated_df), " peaks -> ", arguments$output)
}

main()
