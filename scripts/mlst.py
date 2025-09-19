from mlst_f import write_to_report


with open(snakemake.log[0], "w") as log_file:

    def log(message: str) -> None:
        log_file.write(f"[mlst.py] {message}\n")
        log_file.flush()

    report = snakemake.input[0]
    output = snakemake.output[0]

    log(f"Starting MLST parsing for report={report} -> output={output}")

    write_to_report(report, output, log)

    log(f"Finished MLST parsing for {report}")
