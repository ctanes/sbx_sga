import pandas as pd
from pathlib import Path
import os

def summarize_file(amr_file):
    """Return a dataframe of SNPs annotated with the sample name."""
    if os.path.getsize(amr_file) > 0:
        sample = Path(amr_file).parent.name
        df = pd.read_csv(amr_file, sep="\t")
        df.insert(0, "Sample", sample)
    else:
        df = pd.DataFrame()
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


