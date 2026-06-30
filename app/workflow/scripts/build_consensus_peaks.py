"""Build a consensus peak BED for a single mark.

A peak is kept if it is observed in at least `min_samples` of the mark's
samples: we count, per interval, how many samples overlap it (bedtools
multiinter), keep intervals seen in >= min_samples samples, and merge those
into the final consensus.
"""
import os
import argparse

import pybedtools


def build_consensus(peak_beds, output_path, min_samples):
    """Write a consensus peak BED of intervals observed in >= min_samples samples."""
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # a peak can't be observed in more samples than exist, so cap the threshold
    threshold = min(min_samples, len(peak_beds))

    # multiinter requires sorted inputs; .sort() handles that for us
    sorted_beds = [pybedtools.BedTool(peak_bed).sort() for peak_bed in peak_beds]

    if len(sorted_beds) == 1:
        # only one sample for this mark: consensus is just its merged peaks
        consensus = sorted_beds[0].merge()
    else:
        # column 4 of multiinter output = number of samples overlapping interval
        counts = pybedtools.BedTool().multi_intersect(
            i=[sample_bed.fn for sample_bed in sorted_beds]
        )
        consensus = counts.filter(
            lambda interval: int(interval[3]) >= threshold
        ).merge()

    return consensus.saveas(output_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--beds",
        nargs="+",
        required=True,
        help="Per-sample GoPeaks .bed files belonging to one mark.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the consensus .bed file.",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        required=True,
        help="Minimum number of samples a peak must be observed in to be kept.",
    )
    args = parser.parse_args()

    consensus = build_consensus(args.beds, args.output, args.min_samples)
    print(
        f"samples={len(args.beds)} "
        f"min_samples={min(args.min_samples, len(args.beds))} "
        f"consensus_peaks={consensus.count()}"
    )


if __name__ == "__main__":
    main()
