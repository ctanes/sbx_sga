import os
from functools import reduce
from typing import TextIO

import pandas as pd


def log(message: str, log_file: TextIO) -> None:
    log_file.write(f"[summarize_all.py] {message}\n")
    if hasattr(log_file, "flush"):
        log_file.flush()


def process_filelines(fp, columns, log_file: TextIO):
    log(f"Reading file {fp} with columns {columns}", log_file)
    df = pd.read_csv(fp, sep="\t")
    columns.insert(0, "SampleID")
    log(f"Reordered columns for {fp}: {columns}", log_file)
    return df[columns]


def summarize_all(input_files, output, tools, log_file: TextIO):
    log(
        f"Starting summarize_all for {len(input_files)} files into {output}. "
        f"Available tools: {list(tools.keys())}",
        log_file,
    )
    if not input_files:
        # Handle empty input_files by writing an empty CSV with just the header
        empty_df = pd.DataFrame(columns=["SampleID"])
        empty_df.to_csv(output, index=False, sep="\t")
        log("No input files provided; wrote empty summary", log_file)
        return

    master_list = []
    for fp in input_files:
        tool = os.path.splitext(os.path.basename(fp))[0]
        columns = tools[tool]
        log(f"Processing tool {tool} for file {fp}", log_file)
        df = process_filelines(fp, columns, log_file)
        master_list.append(df)

    final_df = reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"),
        master_list,
    )
    final_df.to_csv(output, index=False, sep="\t")
    log(f"Finished writing combined summary to {output}", log_file)


try:
    log_path = snakemake.log[0]
except (AttributeError, KeyError, IndexError):
    log_path = None

log_handle = open(log_path, "w") if log_path else open(os.devnull, "w")

with log_handle as log_file:
    log("Invoked via Snakemake", log_file)
    summarize_all(snakemake.input, snakemake.output[0], snakemake.params.tools, log_file)
    log("Completed summarize_all Snakemake execution", log_file)
