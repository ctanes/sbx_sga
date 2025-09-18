import sys

from scripts.mash_f import (
    open_report,
    parse_report,
    contamination_call,
    write_report,
)


def log(message: str) -> None:
    sys.stderr.write(f"[mash.py] {message}\n")


sorted_report = snakemake.input[0]
output = snakemake.output[0]

log(f"Starting Mash summary for {sorted_report} -> {output}")

sample, filelines = open_report(sorted_report)
log(f"Loaded {len(filelines)} candidate lines for sample {sample}")

if len(filelines) == 0:
    empty_dict = {"": ""}
    log("No lines found in report; writing empty summary")
    write_report(output, sample, empty_dict)
else:
    parsed_report = parse_report(filelines)
    log(f"Parsed {len(parsed_report)} contamination entries")
    mash_dict = contamination_call(parsed_report)
    log(f"Contamination call result: {mash_dict}")
    write_report(output, sample, mash_dict)

log(f"Finished Mash summary for sample {sample}")
