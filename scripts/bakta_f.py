import sys
import os
import csv


def parse_file(filelines):
    parsed_dict = {}
    if len(filelines) != 0:
        for line in filelines:
            line = line.rstrip().split(":")
            try:
                key = line[0]
                value = line[1].strip()
            except:
                continue
            parsed_dict[key] = value
    return parsed_dict


def test_parse():
    lines = ["Annotation:\n", "test: 9\n", "a: 10\n", "\n", "d: 1343\n"]
    assert parse_file(lines) == {"Annotation": "", "test": "9", "a": "10", "d": "1343"}


def filter_keys(parsed_dict):
    return {key: value for key, value in parsed_dict.items() if value != ""}


def test_filter():
    test = {"Annotation": "", "test": "9", "a": "10", "d": "1343"}
    assert filter_keys(test) == {"test": "9", "a": "10", "d": "1343"}


def write_to_report(report_fp, output_fp):
    with open(report_fp, "r") as f_in:
        lines = f_in.readlines()

    parsed_dict = parse_file(lines)
    filtered = filter_keys(parsed_dict)

    with open(output_fp, "w") as op:
        writer = csv.writer(op, delimiter="\t")
        writer.writerow(filtered.keys())
        writer.writerow(filtered.values())
