from shovill_f import write_shovill_stats

genome = snakemake.input[0]
output = snakemake.output[0]
write_shovill_stats(genome, output)
