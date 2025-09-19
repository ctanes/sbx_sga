import csv
from io import StringIO
from typing import Callable


def _noop_log(message: str) -> None:
    """Default logger used when no log callable is supplied."""
    pass




def parse_file(filelines, log: Callable[[str], None] = _noop_log):
    log(f"[bakta_f] Parsing {len(filelines)} raw lines from Bakta report")
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
            log(f"[bakta_f] Captured entry key={key!r} value={value!r}")
    return parsed_dict


def test_parse():
    lines = ["Annotation:\n", "test: 9\n", "a: 10\n", "\n", "d: 1343\n"]
    buffer = StringIO()

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    assert parse_file(lines, test_log) == {
        "Annotation": "",
        "test": "9",
        "a": "10",
        "d": "1343",
    }


def filter_keys(parsed_dict, log: Callable[[str], None] = _noop_log):
    filtered = {key: value for key, value in parsed_dict.items() if value != ""}
    log(
        "[bakta_f] Filtered parsed entries: "
        f"kept {len(filtered)} of {len(parsed_dict)} keys with non-empty values"
    )
    return filtered


def test_filter():
    test = {"Annotation": "", "test": "9", "a": "10", "d": "1343"}
    buffer = StringIO()

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    assert filter_keys(test, test_log) == {"test": "9", "a": "10", "d": "1343"}


def write_to_report(
    report_fp, output_fp, log: Callable[[str], None] = _noop_log
):
    log(f"[bakta_f] Opening Bakta report at {report_fp}")
    with open(report_fp, "r") as f_in:
        lines = f_in.readlines()

    parsed_dict = parse_file(lines, log)
    filtered = filter_keys(parsed_dict, log)

    with open(output_fp, "w") as op:
        writer = csv.writer(op, delimiter="\t")
        writer.writerow(filtered.keys())
        writer.writerow(filtered.values())
    log(f"[bakta_f] Wrote Bakta summary with {len(filtered)} keys to {output_fp}")
