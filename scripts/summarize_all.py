import os
from functools import reduce
from typing import Callable

import pandas as pd


LogFunc = Callable[[str], None]


def process_filelines(fp, columns, log: LogFunc):
    log(f"[summarize_all] Reading file {fp} with columns {columns}")
    df = pd.read_csv(fp, sep="\t")
    columns.insert(0, "SampleID")
    log(f"[summarize_all] Reordered columns for {fp}: {columns}")
    return df[columns]


def summarize_all(input_files, output, tools, log: LogFunc):
    log(
        "[summarize_all] Starting summarize_all for "
        f"{len(input_files)} files into {output}. "
        f"Available tools: {list(tools.keys())}"
    )
    if not input_files:
        # Handle empty input_files by writing an empty CSV with just the header
        empty_df = pd.DataFrame(columns=["SampleID"])
        empty_df.to_csv(output, index=False, sep="\t")
        log("[summarize_all] No input files provided; wrote empty summary")
        return

    master_list = []
    for fp in input_files:
        tool = os.path.splitext(os.path.basename(fp))[0]
        columns = tools[tool]
        log(f"[summarize_all] Processing tool {tool} for file {fp}")
        df = process_filelines(fp, columns, log)
        master_list.append(df)

    final_df = reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"),
        master_list,
    )
    final_df.to_csv(output, index=False, sep="\t")
    log(f"[summarize_all] Finished writing combined summary to {output}")


with open(snakemake.log[0], "w") as log_file:
    def log(message: str) -> None:
        log_file.write(f"[summarize_all.py] {message}\n")
        log_file.flush()


    log("Invoked via Snakemake")
    summarize_all(snakemake.input, snakemake.output[0], snakemake.params.tools, log)
    log("Completed summarize_all Snakemake execution")
