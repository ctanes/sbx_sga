import pandas as pd
from functools import reduce
import sys


def process_filelines(fp, tool, master_list):
    if tool == "mlst":
        f_obj = open(fp, "r")
        filelines = f_obj.readlines()
        mlst_list = []
        for line in filelines:
            line = line.rstrip().split("\t")
            if len(line) != 10:
                sample = line[0]
                final_line = [sample, "", "", ""]
            else:
                sample = line[0]
                schema = line[1]
                st = line[2]
                alleles = line[3:]
                alleles_joined = " ".join(alleles)
                final_line = [sample, schema, st, alleles_joined]
            mlst_list.append(final_line)
        df = pd.DataFrame(mlst_list)
        df.columns = ["Sample", "Schema", "ST", "Alleles"]
        df["Sample"] = df["Sample"].apply(lambda x: x.rstrip(".fa"))
        master_list.append(df)

    else:
        df = pd.read_csv(fp, sep="\t")
        master_list.append(df)


def summarize_all(input_files, output):
    if not input_files:
        # Handle empty input_files by writing an empty CSV with just the header
        empty_df = pd.DataFrame(columns=["Sample"])
        empty_df.to_csv(output, index=False, sep="\t")
        return
    master_list = []
    for fp in input_files:
        tool = fp.split("/")[-1].split(".report")[0]
        if tool == "amr":
            continue
        else:
            process_filelines(fp, tool, master_list)

    final_df = reduce(
        lambda left, right: pd.merge(left, right, on="Sample", how="outer"), master_list
    )
    final_df.to_csv(output, index=False, sep="\t")



summarize_all(snakemake.input, snakemake.output[0])
