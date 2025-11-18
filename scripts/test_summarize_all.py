import pandas as pd

from scripts.map import antimicrobial, assembly_qc, contaminant, taxonomic_assignment
from scripts.summarize_all import summarize_outputs


def test_summarize_outputs_merges_categories():
    parsed_outputs = {
        "shovill": pd.DataFrame(
            [
                {
                    "SampleID": "S1",
                    "Total_contigs": 10,
                    "Min_coverage": 5,
                    "Max_coverage": 40,
                    "Total_length": 50000,
                    "Average_coverage": 25,
                },
                {
                    "SampleID": "S2",
                    "Total_contigs": 20,
                    "Min_coverage": 4,
                    "Max_coverage": 35,
                    "Total_length": 75000,
                    "Average_coverage": 30,
                },
            ]
        ),
        "bakta": pd.DataFrame(
            [
                {"SampleID": "S1", "GC": 51.1, "N50": 1500, "CDSs": "2000"},
                {"SampleID": "S3", "GC": 47.0, "N50": 900, "CDSs": "1800"},
            ]
        ),
        "checkm": pd.DataFrame(
            [
                {"SampleID": "S2", "Completeness": 98.5, "Contamination": 1.5},
            ]
        ),
        "mash": pd.DataFrame(
            [
                {
                    "SampleID": "S1",
                    "hits_per_thousand": "911/1000",
                    "species": "",
                },
                {
                    "SampleID": "S3",
                    "hits_per_thousand": "110/1000",
                    "species": "Bug",
                },
            ]
        ),
        "mlst": pd.DataFrame(
            [
                {
                    "SampleID": "S1",
                    "classification": "schema1_1 10",
                    "allele_assignment": "geneA(1) geneB(2)",
                },
                {
                    "SampleID": "S3",
                    "classification": "schema2_2 20",
                    "allele_assignment": "geneC(3) geneD(4)",
                },
            ]
        ),
        "sylph": pd.DataFrame(
            [
                {
                    "SampleID": "S2",
                    "Contig_name": "schema_st_15",
                },
                {
                    "SampleID": "S3",
                    "Contig_name": "schema_st_25",
                },
            ]
        ),
        "abritamr": pd.DataFrame(
            [
                {
                    "SampleID": "S2",
                    "Contig id": "contig_1",
                    "Gene symbol": "geneA",
                    "Accession of closest sequence": "acc1",
                    "Element type": "element",
                }
            ]
        ),
    }

    summaries = summarize_outputs(
        parsed_outputs,
        assembly_qc_tools=assembly_qc.keys(),
        taxonomic_assignment_tools=taxonomic_assignment.keys(),
        contaminant_tools=contaminant.keys(),
        antimicrobial_tools=antimicrobial.keys(),
    )

    assembly_qc_df = summaries["assembly_qc"].set_index("SampleID").sort_index()
    assert set(assembly_qc_df.index) == {"S1", "S2", "S3"}
    assert assembly_qc_df.loc["S1", "contig_count"] == 10
    assert assembly_qc_df.loc["S3", "gc_content"] == 47.0
    assert assembly_qc_df.loc["S2", "completeness"] == 98.5

    tax_assignment_df = summaries["taxonomic_assignment"].set_index("SampleID")
    assert set(tax_assignment_df.columns) == {
        "mash_contamination",
        "mash_contaminated_spp",
        "st_schema",
        "st",
        "allele_assignment",
    }
    assert tax_assignment_df.loc["S1", "st_schema"] == "schema1"
    assert tax_assignment_df.loc["S3", "mash_contamination"] == 0.2

    antimicrobial_df = summaries["antimicrobial"]
    assert antimicrobial_df.shape[0] == 1
    assert antimicrobial_df.iloc[0]["contig_id"] == "contig_1"


def test_summarize_outputs_handles_missing_tools():
    summaries = summarize_outputs(
        parsed_outputs={},
        assembly_qc_tools=assembly_qc.keys(),
        taxonomic_assignment_tools=taxonomic_assignment.keys(),
        contaminant_tools=contaminant.keys(),
        antimicrobial_tools=antimicrobial.keys(),
    )

    assert summaries["assembly_qc"].columns.tolist() == ["SampleID"]
    assert summaries["taxonomic_assignment"].columns.tolist() == ["SampleID"]
    assert summaries["antimicrobial"].columns.tolist() == ["SampleID"]
