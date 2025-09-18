import os

from scripts.mash_f import (
    open_report,
    parse_report,
    contamination_call,
    write_report,
)


try:
    log_path = snakemake.log[0]
except (AttributeError, KeyError, IndexError):
    log_path = None

log_handle = open(log_path, "w") if log_path else open(os.devnull, "w")

with log_handle as log_file:
    def log(message: str) -> None:
        log_file.write(f"[mash.py] {message}\n")
        log_file.flush()


    sorted_report = snakemake.input[0]
    output = snakemake.output[0]

    log(f"Starting Mash summary for {sorted_report} -> {output}")

    sample, filelines = open_report(sorted_report, log_file)
    log(f"Loaded {len(filelines)} candidate lines for sample {sample}")

    if len(filelines) == 0:
        empty_dict = {"": ""}
        log("No lines found in report; writing empty summary")
        write_report(output, sample, empty_dict, log_file)
    else:
        parsed_report = parse_report(filelines, log_file)
        log(f"Parsed {len(parsed_report)} contamination entries")
        mash_dict = contamination_call(parsed_report, log_file)
        log(f"Contamination call result: {mash_dict}")
        write_report(output, sample, mash_dict, log_file)

    log(f"Finished Mash summary for sample {sample}")
