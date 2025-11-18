import sys
import pandas as pd
import traceback
from functools import reduce
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.map import (
    antimicrobial,
    assembly_qc,
    reduce_dataframe,
    taxonomic_assignment,
)
from scripts.parse import (
    parse_all_outputs,
    parse_tsv,
    parse_bakta_txt,
    parse_mash_winning_sorted_tab,
    parse_fasta,
)


def _merge_tool_outputs(
    parsed_outputs: dict[str, pd.DataFrame], tools: Iterable[str]
) -> pd.DataFrame:
    dfs = [
        reduce_dataframe(parsed_outputs[tool], tool)
        for tool in tools
        if tool in parsed_outputs
    ]
    if not dfs:
        return pd.DataFrame(columns=["SampleID"])

    return reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"),
        dfs,
    )


def summarize_outputs(
    parsed_outputs: dict[str, pd.DataFrame],
    assembly_qc_tools: Iterable[str],
    taxonomic_assignment_tools: Iterable[str],
    antimicrobial_tools: Iterable[str],
) -> dict[str, pd.DataFrame]:
    return {
        "assembly_qc": _merge_tool_outputs(parsed_outputs, assembly_qc_tools),
        "taxonomic_assignment": _merge_tool_outputs(
            parsed_outputs, taxonomic_assignment_tools
        ),
        "antimicrobial": _merge_tool_outputs(parsed_outputs, antimicrobial_tools),
    }


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log:
        try:
            log.write("Starting summary script\n")

            parsers = {
                "abritamr": parse_tsv,
                "bakta": parse_bakta_txt,
                "checkm": parse_tsv,
                "mash": parse_mash_winning_sorted_tab,
                "mlst": parse_tsv,
                "shovill": parse_fasta,
                "sylph": parse_tsv,
            }

            outputs: dict[str, list[Path]] = {
                "abritamr": [Path(fp) for fp in snakemake.input.abritamr],  # type: ignore
                "bakta": [Path(fp) for fp in snakemake.input.bakta],  # type: ignore
                "checkm": [Path(fp) for fp in snakemake.input.checkm],  # type: ignore
                "mash": [Path(fp) for fp in snakemake.input.mash],  # type: ignore
                "mlst": [Path(fp) for fp in snakemake.input.mlst],  # type: ignore
                "shovill": [Path(fp) for fp in snakemake.input.shovill],  # type: ignore
                "sylph": [Path(fp) for fp in snakemake.input.sylph],  # type: ignore
            }

            tool_reports = {Path(fp).stem: Path(fp) for fp in snakemake.output.tool_reports}  # type: ignore

            assembly_qcs = snakemake.output.assembly_qcs  # type: ignore
            taxonomic_assignments = snakemake.output.taxonomic_assignments  # type: ignore
            antimicrobials = snakemake.output.antimicrobials  # type: ignore

            mash_identity = snakemake.params.mash_identity  # type: ignore
            mash_hits = snakemake.params.mash_hits  # type: ignore
            mash_median_multiplicity_factor = snakemake.params.mash_median_multiplicity_factor  # type: ignore

            parser_kwargs = {
                "mash": {
                    "identity": mash_identity,
                    "hits": mash_hits,
                    "median_multiplicity_factor": mash_median_multiplicity_factor,
                }
            }

            # Parse outputs
            log.write("Parsing tool outputs\n")
            parsed_outputs = parse_all_outputs(outputs, parsers, parser_kwargs)

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
            summary_tables = summarize_outputs(
                parsed_outputs,
                assembly_qc_tools=assembly_qc.keys(),
                taxonomic_assignment_tools=taxonomic_assignment.keys(),
                antimicrobial_tools=antimicrobial.keys(),
            )

            summary_tables["assembly_qc"].to_csv(assembly_qcs, sep="\t", index=False)
            summary_tables["taxonomic_assignment"].to_csv(
                taxonomic_assignments, sep="\t", index=False
            )
            summary_tables["antimicrobial"].to_csv(
                antimicrobials, sep="\t", index=False
            )
            log.write("Finished writing final summaries\n")
        except Exception as error:
            log.write(f"Encountered error: {error}")
            log.write(traceback.format_exc())
            raise
