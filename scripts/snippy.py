import pandas as pd
from pathlib import Path


def summarize_file(snp_file):
    """Return a dataframe of SNPs annotated with the sample name."""
    sample = Path(snp_file).parent.name
    df = pd.read_csv(snp_file, sep="\t")
    df.insert(0, "Sample", sample)
    return df


def summarize_all(files, output):
    frames = [summarize_file(f) for f in files]
    if frames:
        df = pd.concat(frames)
    else:
        df = pd.DataFrame()
    df.to_csv(output, sep="\t", index=False)
    return output

summarize_all(snakemake.input, snakemake.output[0])
