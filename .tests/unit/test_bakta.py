from bakta_f import parse_file, filter_keys, write_to_report


def test_parse_file():
    lines = [
        "Annotation:\n",
        "test: 1\n",
        "foo: 2\n",
    ]
    parsed = parse_file(lines)
    assert parsed == {"Annotation": "", "test": "1", "foo": "2"}


def test_filter_keys():
    parsed = {"Annotation": "", "test": "1", "foo": "2"}
    filtered = filter_keys(parsed)
    assert filtered == {"test": "1", "foo": "2"}


def test_write_to_report(tmp_path):
    report_fp = tmp_path / "report.txt"
    report_fp.write_text("Annotation:\nfoo: 2\nbar: 3\n")

    output_fp = tmp_path / "out.tsv"
    write_to_report(report_fp, output_fp)

    assert output_fp.exists()
    content = output_fp.read_text().splitlines()
    assert content[0] == "foo\tbar"
    assert content[1] == "2\t3"

