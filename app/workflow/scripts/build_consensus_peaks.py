"""Build a consensus peak BED for a single mark.

Run via Snakemake's `script:` directive, so the `snakemake` object is injected
automatically (snakemake.input / output / params / log / wildcards).

A peak is kept if it is observed in at least `min_samples` of the mark's
samples: we count, per interval, how many samples overlap it (bedtools
multiinter), keep intervals seen in >= min_samples samples, and merge those
into the final consensus.
"""
from pathlib import Path
import pybedtools

beds = list(snakemake.input.beds)
out_path = snakemake.output.bed

# a peak can't be observed in more samples than exist, so cap the threshold
min_samples = min(int(snakemake.params.min_samples), len(beds))

Path(out_path).parent.mkdir(parents=True, exist_ok=True)

# multiinter requires sorted inputs; .sort() handles that for us
sorted_beds = [pybedtools.BedTool(bed).sort() for bed in beds]

if len(sorted_beds) == 1:
    # only one sample for this mark: consensus is just its merged peaks
    consensus = sorted_beds[0].merge()
else:
    # column 4 of multiinter output = number of samples overlapping the interval
    counts = pybedtools.BedTool().multi_intersect(i=[b.fn for b in sorted_beds])
    consensus = counts.filter(lambda iv: int(iv[3]) >= min_samples).merge()

consensus.saveas(out_path)

with open(snakemake.log[0], "w") as log:
    log.write(
        f"mark={snakemake.wildcards.mark} "
        f"samples={len(beds)} "
        f"min_samples={min_samples} "
        f"consensus_peaks={pybedtools.BedTool(out_path).count()}\n"
    )
