import pandas as pd
from pathlib import Path
from typing import Callable


def _parse_sample_name(fp: Path) -> str:
    return fp.parent.parent.name


def parse_tsv(fp: Path) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df.insert(0, "SampleID", _parse_sample_name(fp))
    return df


def parse_bakta_txt(fp: Path) -> pd.DataFrame:
    df = pd.DataFrame()
    with open(fp) as f:
        for line in f:
            if len(line.split(":")) != 2:
                continue
            key, value = line.split(":")
            df[key.strip()] = value.strip()

    df.insert(0, "SampleID", _parse_sample_name(fp))
    return df


def parse_mash_winning_sorted_tab(fp: Path) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t", header=None)
    # TODO
    df.insert(0, "SampleID", _parse_sample_name(fp))
    return df


def _parse_header(h: str) -> dict[str, str]:
    hl = h.split()[1:]  # gets rid of the contig name
    header_dict = dict(item.split("=") for item in hl)
    return header_dict


def parse_fasta(fp: Path) -> pd.DataFrame:
    headers = []
    with open(fp) as f:
        for l in f:
            if l.startswith(">"):
                headers.append(l.strip())

    hs = [_parse_header(h) for h in headers]
    total_contigs = len(hs)
    min_cov = min(float(h["cov"]) for h in hs)
    max_cov = max(float(h["cov"]) for h in hs)
    total_length = sum(int(h["len"]) for h in hs)
    avg_cov = sum(int(h["len"]) * float(h["cov"]) for h in hs) / total_length
    rounded_cov = round(avg_cov, 2)

    return pd.DataFrame(
        {
            "SampleID": [_parse_sample_name(fp)],
            "Total_contigs": [total_contigs],
            "Min_coverage": [min_cov],
            "Max_coverage": [max_cov],
            "Total_length": [total_length],
            "Average_coverage": [rounded_cov],
        }
    )


def parse_all_outputs(
    outputs: dict[str, list[Path]], parsers: dict[str, Callable]
) -> dict[str, pd.DataFrame]:
    parsed_outputs = {}
    for tool, fps in outputs.items():
        dfs = [parsers[tool](fp) for fp in fps]
        combined_df = pd.concat(dfs, ignore_index=True)
        parsed_outputs[tool] = combined_df
    return parsed_outputs
