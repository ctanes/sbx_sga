import sys

from mlst_f import write_to_report


def log(message: str) -> None:
    sys.stderr.write(f"[mlst.py] {message}\n")


report = snakemake.input[0]
output = snakemake.output[0]

log(f"Starting MLST parsing for report={report} -> output={output}")

write_to_report(report, output)

log(f"Finished MLST parsing for {report}")
