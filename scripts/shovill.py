import os
import sys


def get_individual_cov(assembly):
    all_ctgs = []
    counter = 1
    for line in assembly:
        line = line.rstrip()
        if line.startswith(">"):
            header_dict = dict(item.split("=") for item in line.split()[1:])
            length = int(header_dict["len"])
            cov = float(header_dict["cov"])
            ctg = counter
            stats = [ctg, cov, length]
            all_ctgs.append(stats)
            counter += 1
    return all_ctgs


def calc_cov_stats(all_ctgs):
    total_ctgs = len(all_ctgs)
    total_length = 0
    total_length_x_cov = 0
    for i in all_ctgs:
        ctg_len = i[-1]
        total_length += ctg_len
        ctg_cov = i[-2]
        len_x_cov = ctg_len * ctg_cov
        total_length_x_cov += len_x_cov
    total_cov = total_length_x_cov / total_length
    rounded_cov = round(total_cov, 2)
    return rounded_cov, total_ctgs


def write_to_report(output, assembly, rounded_cov, total_ctgs):
    sample = os.path.splitext(os.path.basename(assembly))[0]
    with open(f"{output}", "w") as op:
        op.write(f"{sample}\t{rounded_cov}\t{total_ctgs}\n")
    return output



genome = snakemake.input[0]
output = snakemake.output[0]

file_obj = open(genome, "r")
filelines = file_obj.readlines()

if len(filelines) == 0:
    genome_cov = "NA"
    num_ctgs = "NA"
    write_to_report(output, genome, genome_cov, num_ctgs)
else:
    ctg_stats = get_individual_cov(filelines)
    genome_cov, num_ctgs = calc_cov_stats(ctg_stats)
    write_to_report(output, genome, genome_cov, num_ctgs)
