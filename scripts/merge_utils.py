import pandas as pd
from typing import Iterable


def merge_dataframes_on_sampleid(dfs: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Merge dataframes on SampleID while coalescing overlapping columns.

    When multiple tools map to the same target fields (e.g., ``classification``),
    pandas' default merge would create suffixed columns (``_x``/``_y``). This
    helper performs an outer merge and then combines overlapping columns so the
    resulting dataframe retains a single shared column for each field.
    """

    dfs = list(dfs)
    if not dfs:
        return pd.DataFrame(columns=["SampleID"])

    merged = dfs[0]
    for df in dfs[1:]:
        overlapping_cols = [
            col for col in merged.columns if col in df.columns and col != "SampleID"
        ]
        merged = pd.merge(merged, df, on="SampleID", how="outer", suffixes=("", "_dup"))

        for col in overlapping_cols:
            dup_col = f"{col}_dup"
            merged[col] = merged[col].combine_first(merged[dup_col])
            merged = merged.drop(columns=[dup_col])

        dup_cols = [col for col in merged.columns if col.endswith("_dup")]
        if dup_cols:
            merged = merged.drop(columns=dup_cols)

    return merged
