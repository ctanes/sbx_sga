import pandas as pd

from scripts.map import virus as virus_tools
from scripts.summarize_virus import summarize_virus_outputs


def test_summarize_virus_outputs_merges_tools():
    parsed_outputs = {
        "genomad_plasmid_summary": pd.DataFrame(
            [
                {"SampleID": "S1", "seq_name": "plasmid_contig_1", "length": 1000},
                {"SampleID": "S2", "seq_name": "plasmid_contig_2", "length": 2000},
            ]
        ),
        "genomad_virus_summary": pd.DataFrame(
            [
                {"SampleID": "S1", "seq_name": "virus_contig_1", "score": 0.95},
            ]
        ),
        "genomad_plasmid_genes": pd.DataFrame(
            [
                {"SampleID": "S1", "gene_name": "plasmid_gene_1", "function": "foo"},
            ]
        ),
        "genomad_virus_genes": pd.DataFrame(
            [
                {"SampleID": "S2", "contig_id": "virus_gene_contig", "function": "bar"},
            ]
        ),
    }

    summary_df = (
        summarize_virus_outputs(parsed_outputs, virus_tool_names=virus_tools.keys())
        .set_index("SampleID")
        .sort_index()
    )

    assert set(summary_df.index) == {"S1", "S2"}
    assert summary_df.loc["S1", "contig_id_x"] == "plasmid_contig_1"
    assert summary_df.loc["S1", "contig_id_y"] == "virus_contig_1"
    assert summary_df.loc["S1", "gene_x"] == "plasmid_gene_1"
    assert pd.isna(summary_df.loc["S1", "gene_y"])
    assert summary_df.loc["S2", "gene_y"] == "virus_gene_contig"


def test_summarize_virus_outputs_handles_missing_tools():
    summary_df = summarize_virus_outputs(
        parsed_outputs={}, virus_tool_names=virus_tools.keys()
    )

    assert summary_df.columns.tolist() == ["SampleID"]
