import os

from concat_files_f import summarize_all


try:
    log_path = snakemake.log[0]
except (AttributeError, KeyError, IndexError):
    log_path = None

log_handle = open(log_path, "w") if log_path else open(os.devnull, "w")

with log_handle as log_file:
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
        log_file,
    )

    log(f"Finished summarizing files into {output_path}")
