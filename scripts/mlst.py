import sys
from mlst_f import write_to_report

report = snakemake.input[0]
output = snakemake.output[0]

write_to_report(report, output)
