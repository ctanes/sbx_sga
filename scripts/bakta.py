import os

from bakta_f import write_to_report


try:
    log_path = snakemake.log[0]
except (AttributeError, KeyError, IndexError):
    log_path = None

log_handle = open(log_path, "w") if log_path else open(os.devnull, "w")

with log_handle as log_file:
    def log(message: str) -> None:
        log_file.write(f"[bakta.py] {message}\n")
        log_file.flush()


    report = snakemake.input[0]
    output = snakemake.output[0]

    log(f"Starting write_to_report with report={report} -> output={output}")
    write_to_report(report, output, log_file)
    log(f"Finished writing parsed Bakta summary to {output}")
