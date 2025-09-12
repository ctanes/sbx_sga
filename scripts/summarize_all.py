import pandas as pd
from functools import reduce
import os


def process_filelines(fp, columns):
    df = pd.read_csv(fp, sep="\t")
    columns.insert(0, "SampleID")
    return df[columns]


def summarize_all(input_files, output, tools):
    if not input_files:
        # Handle empty input_files by writing an empty CSV with just the header
        empty_df = pd.DataFrame(columns=["SampleID"])
        empty_df.to_csv(output, index=False, sep="\t")
        return

    master_list = []
    for fp in input_files:
        tool = os.path.splitext(os.path.basename(fp))[0]
        columns = tools[tool]
        df = process_filelines(fp, columns)
        master_list.append(df)

    final_df = reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"),
        master_list,
    )
    final_df.to_csv(output, index=False, sep="\t")


summarize_all(snakemake.input, snakemake.output[0], snakemake.params.tools)
