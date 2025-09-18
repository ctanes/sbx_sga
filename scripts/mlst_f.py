import csv
from io import StringIO
from typing import TextIO


def log(message: str, log_file: TextIO) -> None:
    log_file.write(f"[mlst_f] {message}\n")
    if hasattr(log_file, "flush"):
        log_file.flush()


def smush_column(line, log_file: TextIO):
    line_parsed = []
    log(f"Processing MLST line: {line}", log_file)
    if line:
        schema = line[1]
        st = line[2]
        alleles = line[3:]
        alleles_joined = " ".join(alleles)
        line_parsed = [schema, st, alleles_joined]
        log(
            "Extracted schema={schema}, st={st}, allele_count={count}".format(
                schema=schema, st=st, count=len(alleles)
            ),
            log_file,
        )
    return line_parsed


def test_smush_column():
    mlst_row = [
        "marc.bacteremia.1020.a.fa",
        "saureus",
        "30",
        "arcC(2)",
        "aroE(2)",
        "glpF(2)",
        "yqiL(2)",
    ]
    buffer = StringIO()
    assert smush_column(mlst_row, buffer) == [
        "saureus",
        "30",
        "arcC(2) aroE(2) glpF(2) yqiL(2)",
    ]


def write_to_report(report, output, log_file: TextIO):
    log(f"Writing MLST summary from {report} to {output}", log_file)
    with open(report, "r") as f_in:
        reader = csv.reader(f_in, delimiter="\t")
        with open(output, "w") as op:
            writer = csv.writer(op, delimiter="\t")
            writer.writerow(["Schema", "ST", "Alleles"])
            for row in reader:
                smushed_column = smush_column(row, log_file)
                writer.writerow(smushed_column)
    log(f"Completed MLST summary for {report}", log_file)
