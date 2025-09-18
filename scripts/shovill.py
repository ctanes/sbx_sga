import sys

from shovill_f import write_shovill_stats


def log(message: str) -> None:
    sys.stderr.write(f"[shovill.py] {message}\n")


genome = snakemake.input[0]
output = snakemake.output[0]

log(f"Starting Shovill stats calculation for genome={genome} -> output={output}")
write_shovill_stats(genome, output)
log(f"Finished Shovill stats calculation for {genome}")
