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
