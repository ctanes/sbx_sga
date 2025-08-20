import sys
import os
from pathlib import Path


# Parsing the tsv file
def parse_file(filepath):
    with open(filepath, "r") as file_obj:
        filelines = file_obj.readlines()
        if len(filelines) >= 2:
            keys = filelines[0].strip().split("\t")
            hits = []
            for line in filelines[1:]:
                values = line.strip().split("\t")
                hit_dict = dict(zip(keys, values))
                hits.append(hit_dict)
            return hits
        else:
            return [
                {
                    "Taxonomic_abundance": "NA",
                    "Contig_name": "NA",
                }
            ]


# Fetching Sample name, taxonomic abundance, and contig name from parsed data
def get_stats(data_dict, report):
    sample_name = os.path.basename(report).split(".tsv")[0]
    taxo_abundance = data_dict.get("Taxonomic_abundance", "NA")
    contig = data_dict.get("Contig_name", "NA")
    return sample_name, taxo_abundance, contig


# Writing it to the snakemake output
def write_report(output, hits, report):
    with open(output, "w") as op:
        for hit in hits:
            sample_name, taxo_abundance, contig = get_stats(hit, report)
            op.write(f"{sample_name}\t{taxo_abundance}\t{contig}\n")
    return output




sys.stderr.write("Starting Sylph processing...\n")

report = snakemake.input[0]
output = snakemake.output[0]

sys.stderr.write(f"Processing report: {report}\n")
sylph_stats = parse_file(report)
sys.stderr.write(f"Parsed data: {sylph_stats}\n")
sys.stderr.write("Obtained stats!\n")
write_report(output, sylph_stats, report)
sys.stderr.write(f"Report written to: {output}\n")
