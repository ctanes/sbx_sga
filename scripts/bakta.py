import sys

from bakta_f import write_to_report


def log(message: str) -> None:
    sys.stderr.write(f"[bakta.py] {message}\n")


report = snakemake.input[0]
output = snakemake.output[0]

log(f"Starting write_to_report with report={report} -> output={output}")
write_to_report(report, output)
log(f"Finished writing parsed Bakta summary to {output}")
