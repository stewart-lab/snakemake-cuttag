"""Build a raw count matrix for a single mark.

Counts, for each consensus peak of the mark, how many reads from each sample's
filtered .bam overlap that peak (via `bedtools multicov`). The result is a
peaks-by-samples matrix of raw integer counts, which feeds the downstream
sample-dropping, low-count filtering, PCA and DESeq2 steps.

Rows are consensus peaks (identified as "chrom:start-end"), columns are samples.
The peak order matches the consensus .bed; the sample column order matches the
order in which --bams / --samples were given.
"""
import os
import csv
import argparse

import pybedtools


def build_count_matrix(consensus_bed, bams):
    """Return (peak_ids, counts) for the consensus peaks against each bam.

    peak_ids is a list of "chrom:start-end" strings (one per consensus peak).
    counts is a list of rows, each row holding one integer count per bam, in
    the same order as `bams`.
    """
    consensus = pybedtools.BedTool(consensus_bed)

    # `bedtools multicov` appends one count column per bam, in bam order. The
    # bams must be indexed (.bai); the filtered-reads step produces those.
    covered = consensus.multi_bam_coverage(bams=bams)

    n_bams = len(bams)
    peak_ids, counts = [], []
    for interval in covered:
        peak_ids.append(f"{interval.chrom}:{interval.start}-{interval.end}")
        # the appended counts are the last n_bams fields of the interval
        fields = interval.fields
        counts.append([int(value) for value in fields[-n_bams:]])

    return peak_ids, counts


def write_matrix(output_path, peak_ids, sample_names, counts):
    """Write the peaks-by-samples matrix as CSV with a 'peak' index column."""
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak"] + sample_names)
        for peak_id, row in zip(peak_ids, counts):
            writer.writerow([peak_id] + row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--peaks",
        required=True,
        help="Consensus peak .bed file for one mark.",
    )
    parser.add_argument(
        "--bams",
        nargs="+",
        required=True,
        help="Filtered, indexed .bam files (one per sample of the mark).",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample names matching the order of --bams (used as column headers).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the raw count matrix .csv.",
    )
    args = parser.parse_args()

    if len(args.bams) != len(args.samples):
        parser.error(
            f"--bams ({len(args.bams)}) and --samples ({len(args.samples)}) "
            "must have the same length"
        )

    peak_ids, counts = build_count_matrix(args.peaks, args.bams)
    write_matrix(args.output, peak_ids, args.samples, counts)
    print(
        f"samples={len(args.samples)} "
        f"peaks={len(peak_ids)} "
        f"matrix={args.output}"
    )


if __name__ == "__main__":
    main()
