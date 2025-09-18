import os

from shovill_f import write_shovill_stats


try:
    log_path = snakemake.log[0]
except (AttributeError, KeyError, IndexError):
    log_path = None

log_handle = open(log_path, "w") if log_path else open(os.devnull, "w")

with log_handle as log_file:
    def log(message: str) -> None:
        log_file.write(f"[shovill.py] {message}\n")
        log_file.flush()


    genome = snakemake.input[0]
    output = snakemake.output[0]

    log(f"Starting Shovill stats calculation for genome={genome} -> output={output}")
    write_shovill_stats(genome, output, log_file)
    log(f"Finished Shovill stats calculation for {genome}")
