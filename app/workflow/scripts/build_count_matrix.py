"""Build a raw count matrix for a single mark.

For each consensus peak of the mark, count how many sequenced fragments from
each sample's filtered .bam overlap that peak, using featureCounts in
paired-end fragment-counting mode (`-p --countReadPairs`). CUT&Tag is
paired-end, so counting fragments (read pairs) rather than individual reads is
the biologically meaningful unit and what DESeq2 expects.

The result is a peaks-by-samples matrix of raw integer counts, which feeds the
downstream sample-dropping, low-count filtering, PCA and DESeq2 steps. Rows are
consensus peaks (identified as "chrom:start-end", in BED 0-based coordinates),
columns are samples in the order --bams / --samples were given.
"""
import os
import csv
import argparse
import tempfile
import subprocess


def read_consensus_peaks(consensus_bed):
    """Yield (peak_id, chrom, start, end) for each interval in a BED file.

    peak_id is "chrom:start-end" in the BED's native 0-based half-open
    coordinates; start/end are returned as ints. Header/track lines skipped.
    """
    with open(consensus_bed) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            start, end = int(start), int(end)
            yield f"{chrom}:{start}-{end}", chrom, start, end


def write_saf(peaks, saf_path):
    """Write peaks to a featureCounts SAF annotation file.

    SAF coordinates are 1-based and inclusive, so the BED 0-based half-open
    [start, end) becomes [start + 1, end]. Strand is unused (we count
    unstranded) but the column is required, so we emit "+".
    """
    with open(saf_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
        for peak_id, chrom, start, end in peaks:
            writer.writerow([peak_id, chrom, start + 1, end, "+"])


def run_featurecounts(saf_path, bams, counts_path, threads):
    """Run featureCounts over the SAF peaks for all bams (fragment counting).

    featureCounts' progress/summary goes to stderr, which it shares with this
    process; the Snakemake rule redirects both to the rule log.
    """
    subprocess.run(
        [
            "featureCounts",
            "-F", "SAF",
            "-a", saf_path,
            "-o", counts_path,
            "-p", "--countReadPairs",  # count fragments (read pairs), not reads
            "-T", str(threads),
            *bams,
        ],
        check=True,
    )


def parse_featurecounts(counts_path, n_bams):
    """Parse a featureCounts output table into (peak_ids, counts).

    The table has a leading '#' comment line, then a header row
    (Geneid, Chr, Start, End, Strand, Length, <bam1>, <bam2>, ...) followed by
    one row per feature. The per-bam counts are the last n_bams columns, in the
    order the bams were passed on the command line.
    """
    peak_ids, counts = [], []
    with open(counts_path) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            fields = line.rstrip("\n").split("\t")
            peak_ids.append(fields[0])
            counts.append([int(value) for value in fields[-n_bams:]])
    return peak_ids, counts


def write_matrix(output_path, peak_ids, sample_names, counts):
    """Write the peaks-by-samples matrix as CSV with a 'peak' index column."""
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak"] + sample_names)
        for peak_id, row in zip(peak_ids, counts):
            writer.writerow([peak_id] + row)


def build_count_matrix(consensus_bed, bams, sample_names, output_path, threads):
    """Count fragments per consensus peak per sample and write the CSV matrix."""
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    peaks = list(read_consensus_peaks(consensus_bed))

    # featureCounts errors on an empty annotation; an empty consensus simply
    # yields a header-only matrix (zero peaks), which is a valid degenerate case.
    if not peaks:
        write_matrix(output_path, [], sample_names, [])
        return 0

    with tempfile.TemporaryDirectory() as tmp_dir:
        saf_path = os.path.join(tmp_dir, "peaks.saf")
        counts_path = os.path.join(tmp_dir, "counts.txt")
        write_saf(peaks, saf_path)
        run_featurecounts(saf_path, bams, counts_path, threads)
        peak_ids, counts = parse_featurecounts(counts_path, len(bams))

    write_matrix(output_path, peak_ids, sample_names, counts)
    return len(peak_ids)


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
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for featureCounts.",
    )
    args = parser.parse_args()

    if len(args.bams) != len(args.samples):
        parser.error(
            f"--bams ({len(args.bams)}) and --samples ({len(args.samples)}) "
            "must have the same length"
        )

    n_peaks = build_count_matrix(
        args.peaks, args.bams, args.samples, args.output, args.threads
    )
    print(
        f"samples={len(args.samples)} "
        f"peaks={n_peaks} "
        f"matrix={args.output}"
    )


if __name__ == "__main__":
    main()
