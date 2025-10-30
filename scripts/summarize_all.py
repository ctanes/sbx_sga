import pandas as pd
import traceback
from pathlib import Path
from .tools import (
    BaktaOutput,
    CheckMOutput,
    MashOutput,
    MLSTOutput,
    ShovillOutput,
    SylphOutput,
    ToolOutput,
)


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log:
        try:
            log.write("Starting summary script\n")

            reports: list[tuple[Path, type[ToolOutput]]] = [
                *[(Path(i), BaktaOutput) for i in snakemake.input.bakta],  # type: ignore
                *[(Path(i), CheckMOutput) for i in snakemake.input.checkm],  # type: ignore
                *[(Path(i), MashOutput) for i in snakemake.input.mash],  # type: ignore
                *[(Path(i), MLSTOutput) for i in snakemake.input.mlst],  # type: ignore
                *[(Path(i), ShovillOutput) for i in snakemake.input.shovill],  # type: ignore
                *[(Path(i), SylphOutput) for i in snakemake.input.sylph],  # type: ignore
            ]
            samples = {fp.parent.name for fp, _ in reports}
            tool_outputs_per_sample: dict[str, list[ToolOutput]] = {}
            for sample in samples:
                sample_tool_outputs = []
                for fp, tool_cls in reports:
                    if fp.parent.name != sample:
                        continue
                    output = tool_cls.from_report(fp, log.write)
                    sample_tool_outputs.append(output)
                tool_outputs_per_sample[sample] = sample_tool_outputs

            final_summary: list[dict[str, str]] = []
            for sample, tool_outputs in tool_outputs_per_sample.items():
                entry = {"SampleID": sample}
                for to in tool_outputs:
                    for key in to.KEYS:
                        entry[key] = to.d.get(key, "NA")
                final_summary.append(entry)

            out_fp = snakemake.output[0]  # type: ignore
            df = pd.DataFrame(final_summary)
            df = df.sort_values("SampleID")
            df.to_csv(out_fp, sep="\t", index=False)
            log.write(f"Wrote final summary to {out_fp}\n")
        except Exception as error:
            log.write(f"Encountered error: {error}")
            log.write(traceback.format_exc())
            raise
