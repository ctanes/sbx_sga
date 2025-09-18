from concat_files_f import summarize_all


with open(snakemake.log[0], "w") as log_file:
    def log(message: str) -> None:
        log_file.write(f"[concat_files.py] {message}\n")
        log_file.flush()


    input_files = list(snakemake.input)
    output_path = snakemake.output[0]
    suffix = snakemake.params.suffix
    header = snakemake.params.header

    log(
        "Starting summarize_all with "
        f"{len(input_files)} files, suffix={suffix!r}, header={header} -> {output_path}"
    )

    summarize_all(
        input_files,
        output_path,
        suffix,
        header,
        log,
    )

    log(f"Finished summarizing files into {output_path}")
