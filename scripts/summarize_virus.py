import sys
import pandas as pd
import traceback
from functools import reduce
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.map import reduce_dataframe, virus as virus_tools


def summarize_virus_outputs(
    parsed_outputs: dict[str, pd.DataFrame], virus_tool_names: Iterable[str]
) -> pd.DataFrame:
    dfs = [
        reduce_dataframe(parsed_outputs[tool], tool)
        for tool in virus_tool_names
        if tool in parsed_outputs
    ]
    if not dfs:
        return pd.DataFrame(columns=["SampleID"])

    return reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"),
        dfs,
    )


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log:
        try:
            log.write("Starting summary script\n")
            from scripts.parse import parse_all_outputs, parse_tsv

            parsers = {
                "genomad_plasmid_summary": parse_tsv,
                "genomad_virus_summary": parse_tsv,
                "genomad_plasmid_genes": parse_tsv,
                "genomad_virus_genes": parse_tsv,
            }

            outputs: dict[str, list[Path]] = {
                "genomad_plasmid_summary": [Path(fp) for fp in snakemake.input if "plasmid_summary" in fp],  # type: ignore
                "genomad_virus_summary": [Path(fp) for fp in snakemake.input if "virus_summary" in fp],  # type: ignore
                "genomad_plasmid_genes": [Path(fp) for fp in snakemake.input if "plasmid_genes" in fp],  # type: ignore
                "genomad_virus_genes": [Path(fp) for fp in snakemake.input if "virus_genes" in fp],  # type: ignore
            }

            tool_reports = {Path(fp).stem: Path(fp) for fp in snakemake.output.tool_reports}  # type: ignore

            virus_summary_output = Path(snakemake.output.virus)  # type: ignore

            # Parse outputs
            log.write("Parsing tool outputs\n")
            parsed_outputs = parse_all_outputs(outputs, parsers)

            # Write individual tool reports
            log.write("Writing individual tool reports\n")
            if set(tool_reports.keys()) != set(parsed_outputs.keys()):
                log.write(
                    f"Warning: tool reports keys {list(tool_reports.keys())} do not match parsed outputs keys {list(parsed_outputs.keys())}\n"
                )
            for tool, df in parsed_outputs.items():
                df.to_csv(tool_reports[tool], sep="\t", index=False)

            # Produce final summaries
            log.write("Producing final summaries\n")
            virus_df = summarize_virus_outputs(
                parsed_outputs, virus_tool_names=virus_tools.keys()
            )
            virus_df.to_csv(virus_summary_output, sep="\t", index=False)
            log.write("Finished writing final summaries\n")
        except Exception as error:
            log.write(f"Encountered error: {error}")
            log.write(traceback.format_exc())
            raise
