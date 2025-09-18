import pandas as pd
from functools import reduce
import os
import sys


def log(message: str) -> None:
    sys.stderr.write(f"[summarize_all.py] {message}\n")


def process_filelines(fp, columns):
    log(f"Reading file {fp} with columns {columns}")
    df = pd.read_csv(fp, sep="\t")
    columns.insert(0, "SampleID")
    log(f"Reordered columns for {fp}: {columns}")
    return df[columns]


def summarize_all(input_files, output, tools):
    log(
        f"Starting summarize_all for {len(input_files)} files into {output}. "
        f"Available tools: {list(tools.keys())}"
    )
    if not input_files:
        # Handle empty input_files by writing an empty CSV with just the header
        empty_df = pd.DataFrame(columns=["SampleID"])
        empty_df.to_csv(output, index=False, sep="\t")
        log("No input files provided; wrote empty summary")
        return

    master_list = []
    for fp in input_files:
        tool = os.path.splitext(os.path.basename(fp))[0]
        columns = tools[tool]
        log(f"Processing tool {tool} for file {fp}")
        df = process_filelines(fp, columns)
        master_list.append(df)

    final_df = reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"),
        master_list,
    )
    final_df.to_csv(output, index=False, sep="\t")
    log(f"Finished writing combined summary to {output}")


log("Invoked via Snakemake")
summarize_all(snakemake.input, snakemake.output[0], snakemake.params.tools)
log("Completed summarize_all Snakemake execution")
