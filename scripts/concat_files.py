from concat_files_f import summarize_all

summarize_all(snakemake.input, snakemake.output[0], snakemake.params.suffix, snakemake.params.header)
