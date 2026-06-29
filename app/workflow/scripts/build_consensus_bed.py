import os

def build_consensus_bed(beds: list[str], output_path: str, n_required_observed: int = 2):
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    n_required_observed = min(len(beds), n_required_observed)

    # output path must end in .bed
    if not output_path.endswith(".bed"):
        raise ValueError("Output file must end with .bed")
    
    # check if all bed files exist
    for bed in beds:
        if not os.path.exists(bed):
            raise FileNotFoundError(f"Bed file {bed} does not exist.")

    # calculate intervals in common with 'bedtools multiinter'
    # https://bedtools.readthedocs.io/en/latest/content/tools/multiinter.html
    interval_file_path = output_path.replace(".bed", ".intervals")

    cmd = [
        "bedtools", "multiinter",
        "-i", *beds,
        ">", interval_file_path,
    ]

    cmd_str = " ".join(cmd)
    os.system(cmd_str)

    # get the intervals that are in at least n_required_observed samples
    with open(interval_file_path, "r") as f:
        lines = f.readlines()

    unmerged_consensus_file = output_path.replace(".bed", "_unmerged.bed")
    output = []
    for line in lines:
        spl = line.strip().split("\t")
        interval_count = int(spl[3])

        if interval_count >= n_required_observed:
            bed_line = "\t".join(spl[:3])
            output.append(bed_line)

    with open(unmerged_consensus_file, "w") as f:
        for line in output:
            f.write(line + "\n")

    # use 'bedtools merge' to consolidate the peaks
    # https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
    cmd = [
        "bedtools", "merge",
        "-i", unmerged_consensus_file,
        ">", output_path,
    ]
    cmd_str = " ".join(cmd)
    os.system(cmd_str)