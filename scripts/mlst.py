import traceback

from mlst_f import write_to_report


with open(snakemake.log[0], "w") as log_file:

    def log(message: str) -> None:
        log_file.write(f"[mlst.py] {message}\n")
        log_file.flush()

    report = snakemake.input[0]
    output = snakemake.output[0]

    try:
        log(f"Starting MLST parsing for report={report} -> output={output}")

        write_to_report(report, output, log)

        log(f"Finished MLST parsing for {report}")
    except Exception as error:
        log(f"Encountered error: {error}")
        log(traceback.format_exc())
        raise
