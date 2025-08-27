import sys
from bakta_f import parse_file, get_annotation_stats, write_to_report

report = snakemake.input[0]
output = snakemake.output[0]

filelines = open(report, "r").readlines()
bakta_dict = parse_file(filelines)
write_to_report(output, report, bakta_dict)
