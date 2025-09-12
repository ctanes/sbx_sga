import pytest
pandas = pytest.importorskip("pandas")
from summarize_all import process_filelines, summarize_all


def test_process_filelines(tmp_path):
    input_file = tmp_path / "test_input.txt"
    input_file.write_text("SampleID\tA\tB\nsample1\t1\t2\n")
    df = process_filelines(input_file, ["A", "B"])
    assert list(df.columns) == ["SampleID", "A", "B"]
    assert df.iloc[0].tolist() == ["sample1", 1, 2]


def test_summarize_all(tmp_path):
    file1 = tmp_path / "tool1.tsv"
    file1.write_text("SampleID\tA\ns1\t1\n")
    file2 = tmp_path / "tool2.tsv"
    file2.write_text("SampleID\tB\ns1\t2\n")
    output = tmp_path / "summary.tsv"
    tools = {"tool1": ["A"], "tool2": ["B"]}
    summarize_all([file1, file2], output, tools)
    df = pandas.read_csv(output, sep="\t")
    assert list(df.columns) == ["SampleID", "A", "B"]
    assert df.iloc[0].tolist() == ["s1", 1, 2]

