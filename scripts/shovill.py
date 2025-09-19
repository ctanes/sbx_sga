from shovill_f import write_shovill_stats


with open(snakemake.log[0], "w") as log_file:

    def log(message: str) -> None:
        log_file.write(f"[shovill.py] {message}\n")
        log_file.flush()

    genome = snakemake.input[0]
    output = snakemake.output[0]

    log(f"Starting Shovill stats calculation for genome={genome} -> output={output}")
    write_shovill_stats(genome, output, log)
    log(f"Finished Shovill stats calculation for {genome}")
