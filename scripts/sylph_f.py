import sys
import os
from pathlib import Path

# Parsing the tsv file
def parse_file(filepath):
    with open(filepath, "r") as file_obj:
        filelines = file_obj.readlines()
        if len(filelines) > 2:
            keys = filelines[0].strip().split("\t")
            values = filelines[1].strip().split("\t")
            data_dict = dict(zip(keys, values))
            return data_dict
        else:
            return {
                "Sample_file": filepath.split(".tsv")[0],
                "Taxonomic_abundance": "NA",
                "Contig_name": "NA",
            }


# Fetching Sample name, taxonomic abundance, and contig name from parsed data
def get_stats(data_dict):
    sample_name = data_dict.get("Sample_file", "NA")
    taxo_abundance = data_dict.get("Taxonomic_abundance", "NA")
    contig = data_dict.get("Contig_name", "NA")
    return sample_name, taxo_abundance, contig


# Writing it to the snakemake output
def write_report(output, sample_name, taxo_abundance, contig):
    with open(output, "w") as op:
        sample_name = os.path.basename(sample_name).split("_1.fastq.gz")[0]
        op.write(f"{sample_name}\t{taxo_abundance}\t{contig}\n")
    return output
