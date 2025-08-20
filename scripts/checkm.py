import sys
import os
from pathlib import Path


def parse_file(quality_report):
    sample_path = Path(quality_report)
    sample_name = sample_path.parts[-2]

    with open(quality_report, "r") as file_obj:
        filelines = file_obj.readlines()

        if len(filelines) != 0:
            keys = filelines[0].strip().split("\t")
            values = filelines[1].strip().split("\t")
            data_dict = dict(zip(keys, values))
            return data_dict
        else:
            return {"Name": sample_name, "Completeness": "NA", "Contamination": "NA"}


def get_stats(data_dict):
    completeness = data_dict["Completeness"]
    contamination = data_dict["Contamination"]
    sample = data_dict["Name"]
    return sample, completeness, contamination


def write_report(output, sample, completeness, contamination):
    with open(output, "w") as op:
        op.write(f"{sample}\t{completeness}\t{contamination}\n")
    return output

report = snakemake.input[0]
output = snakemake.output[0]

checkm_stats = parse_file(report)
sample, completeness, contamination = get_stats(checkm_stats)
write_report(output, sample, completeness, contamination)
