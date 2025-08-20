import csv
from pathlib import Path
import sys

def parse_file(file):
    sample = Path(file)
    sample_name = sample.parts[-2]
    with open(file, "r") as file_obj:
        filelines = file_obj.readlines()
    return sample_name, filelines


# Expects a list of files as the argument
def parse_scope(files):
    subclassList = {}

    for file in files:
        sample, filelines = parse_file(file)
        for line in filelines[1:]:
            line = line.rstrip()
            line_list = line.split("\t")
            scope = line_list[8]
            gene = line_list[5]

            if scope == "AMR":
                subclass = line_list[11]
                if subclass not in subclassList:
                    subclassList[subclass] = {}
                if sample not in subclassList[subclass]:
                    subclassList[subclass][sample] = [gene]
                else:
                    subclassList[subclass][sample].append(gene)
            else:
                continue
    return subclassList


def write_long_report(subclassList, output):
    with open(output, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(["Sample", "Gene", "Category"])

        for category, samples in subclassList.items():
            for sample, genes in samples.items():
                for gene in genes:
                    writer.writerow([sample, gene, category])
    return output



input_reports = snakemake.input
output = snakemake.output[0]

subclass_list = parse_scope(input_reports)
write_long_report(subclass_list, output)
