import csv
from io import StringIO
import os
from concat_files_f import write_report, summarize_all


def test_write_report():
    input_io = StringIO("col1\tcol2\nval1\tval2\n")
    reader = csv.reader(input_io, delimiter="\t")
    output_io = StringIO()
    writer = csv.writer(output_io, delimiter="\t")

    # skip header
    next(reader)
    write_report(writer, reader, "sampleA")

    output_io.seek(0)
    lines = output_io.getvalue().strip().splitlines()
    assert lines == ["sampleA\tval1\tval2"]


def test_summarize_all(tmp_path):
    sample1 = tmp_path / "sample1"
    sample1.mkdir()
    report1 = sample1 / "report.tsv"
    report1.write_text("col1\tcol2\n1\t2\n")

    sample2 = tmp_path / "sample2"
    sample2.mkdir()
    report2 = sample2 / "report.tsv"
    report2.write_text("col1\tcol2\n3\t4\n")

    output = tmp_path / "out.tsv"
    summarize_all([report1, report2], output)

    lines = output.read_text().splitlines()
    assert lines[0] == "SampleID\tcol1\tcol2"
    assert lines[1] == "sample1\t1\t2"
    assert lines[2] == "sample2\t3\t4"

