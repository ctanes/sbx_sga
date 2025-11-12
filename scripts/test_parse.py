import pytest
from pathlib import Path
from .parse import (
    parse_tsv,
    parse_bakta_txt,
    parse_mash_winning_sorted_tab,
    parse_fasta,
)


@pytest.fixture
def test_reports_fp():
    return Path(__file__).parent.parent / ".tests/data/test_reports"


def test_parse_mash_marc_3111(test_reports_fp):
    fp = test_reports_fp / "mash/marc.bacteremia.3111_sorted_winning.tab"
    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )
    assert not df.empty
    assert "species" in df.columns
    assert all(df["identity"] >= 0.85)
    assert all(df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= 100)
    assert df["species"].iloc[0] == "Enterobacter cloacae"
    assert df["species"].iloc[1] == "Enterobacter kobei"


def test_parse_mash_s234_ori(test_reports_fp):
    fp = test_reports_fp / "mash/s234.ori.lightblue.b_sorted_winning.tab"
    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )
    assert not df.empty
    assert all(df["identity"] >= 0.85)
    assert all(df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= 100)
    assert df["species"].iloc[0] == "Bacillus cereus"
    assert df["species"].iloc[1] == "Pseudomonas denitrificans"
    assert df["species"].iloc[2] == "Stenotrophomonas maltophilia"
    assert df["species"].iloc[3] == "Stenotrophomonas acidaminiphila"


def test_parse_mash_marc_235(test_reports_fp):
    fp = test_reports_fp / "mash/marc.entero.235_sorted_winning.tab"
    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )
    assert df.empty

    df = parse_mash_winning_sorted_tab(
        fp, identity=0.80, hits=10, median_multiplicity_factor=0.05
    )
    assert not df.empty
    assert all(df["identity"] >= 0.80)
    assert all(df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= 10)
    assert df["species"].iloc[0] == "Serratia marcescens"


def test_fasta(test_reports_fp):
    fp = test_reports_fp / "shovill/dummy.fa"
    df = parse_fasta(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Total_contigs" in df.columns


def test_empty_fasta(test_reports_fp):
    fp = test_reports_fp / "shovill/empty.fa"
    df = parse_fasta(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert df["Total_contigs"].iloc[0] == 0


def test_bakta(test_reports_fp):
    fp = test_reports_fp / "bakta/dummy.txt"
    df = parse_bakta_txt(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "CDSs" in df.columns
    assert df["CDSs"].iloc[0] == "2492"


def test_parse_tsv_sylph_3151(test_reports_fp):
    fp = test_reports_fp / "sylph/marc.bacteremia.3151.tsv"
    df = parse_tsv(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Contig_name" in df.columns
    assert "Actinomyces" in df["Contig_name"].iloc[0]
