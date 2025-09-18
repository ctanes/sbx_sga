import csv
from io import StringIO
from typing import TextIO


def log(message: str, log_file: TextIO) -> None:
    log_file.write(f"[bakta_f] {message}\n")
    if hasattr(log_file, "flush"):
        log_file.flush()


def parse_file(filelines, log_file: TextIO):
    log(f"Parsing {len(filelines)} raw lines from Bakta report", log_file)
    parsed_dict = {}
    if len(filelines) != 0:
        for line in filelines:
            line = line.rstrip().split(":")
            try:
                key = line[0]
                value = line[1].strip()
            except Exception:
                continue
            parsed_dict[key] = value
            log(f"Captured entry key={key!r} value={value!r}", log_file)
    return parsed_dict


def test_parse():
    lines = ["Annotation:\n", "test: 9\n", "a: 10\n", "\n", "d: 1343\n"]
    buffer = StringIO()
    assert parse_file(lines, buffer) == {
        "Annotation": "",
        "test": "9",
        "a": "10",
        "d": "1343",
    }


def filter_keys(parsed_dict, log_file: TextIO):
    filtered = {key: value for key, value in parsed_dict.items() if value != ""}
    log(
        "Filtered parsed entries: "
        f"kept {len(filtered)} of {len(parsed_dict)} keys with non-empty values",
        log_file,
    )
    return filtered


def test_filter():
    test = {"Annotation": "", "test": "9", "a": "10", "d": "1343"}
    buffer = StringIO()
    assert filter_keys(test, buffer) == {"test": "9", "a": "10", "d": "1343"}


def write_to_report(report_fp, output_fp, log_file: TextIO):
    log(f"Opening Bakta report at {report_fp}", log_file)
    with open(report_fp, "r") as f_in:
        lines = f_in.readlines()

    parsed_dict = parse_file(lines, log_file)
    filtered = filter_keys(parsed_dict, log_file)

    with open(output_fp, "w") as op:
        writer = csv.writer(op, delimiter="\t")
        writer.writerow(filtered.keys())
        writer.writerow(filtered.values())
    log(
        f"Wrote Bakta summary with {len(filtered)} keys to {output_fp}",
        log_file,
    )
