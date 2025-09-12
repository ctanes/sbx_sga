from mlst_f import smush_column, write_to_report


def test_smush_column():
    row = ["sample.fa", "saureus", "30", "arcC(2)", "aroE(2)"]
    assert smush_column(row) == ["saureus", "30", "arcC(2) aroE(2)"]


def test_write_to_report(tmp_path):
    report = tmp_path / "mlst.tsv"
    report.write_text("sample.fa\tsaureus\t30\tarcC(2)\taroE(2)\n")
    output = tmp_path / "out.tsv"
    write_to_report(report, output)

    lines = output.read_text().splitlines()
    assert lines[0] == "Schema\tST\tAlleles"
    assert lines[1] == "saureus\t30\tarcC(2) aroE(2)"

