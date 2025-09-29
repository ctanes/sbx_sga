import traceback

from bakta_f import write_to_report


with open(snakemake.log[0], "w") as log_file:

    def log(message: str) -> None:
        log_file.write(f"[bakta.py] {message}\n")
        log_file.flush()

    report = snakemake.input[0]
    output = snakemake.output[0]

    try:
        log(f"Starting write_to_report with report={report} -> output={output}")
        write_to_report(report, output, log)
        log(f"Finished writing parsed Bakta summary to {output}")
    except Exception as error:
        log(f"Encountered error: {error}")
        log(traceback.format_exc())
        raise
