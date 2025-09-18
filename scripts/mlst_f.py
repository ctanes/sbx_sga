import csv
from io import StringIO
from typing import Callable


LogFunc = Callable[[str], None]


def smush_column(line, log: LogFunc):
    line_parsed = []
    log(f"[mlst_f] Processing MLST line: {line}")
    if line:
        schema = line[1]
        st = line[2]
        alleles = line[3:]
        alleles_joined = " ".join(alleles)
        line_parsed = [schema, st, alleles_joined]
        log(
            "[mlst_f] Extracted schema={schema}, st={st}, allele_count={count}".format(
                schema=schema, st=st, count=len(alleles)
            )
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

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    assert smush_column(mlst_row, test_log) == [
        "saureus",
        "30",
        "arcC(2) aroE(2) glpF(2) yqiL(2)",
    ]


def write_to_report(report, output, log: LogFunc):
    log(f"[mlst_f] Writing MLST summary from {report} to {output}")
    with open(report, "r") as f_in:
        reader = csv.reader(f_in, delimiter="\t")
        with open(output, "w") as op:
            writer = csv.writer(op, delimiter="\t")
            writer.writerow(["Schema", "ST", "Alleles"])
            for row in reader:
                smushed_column = smush_column(row, log)
                writer.writerow(smushed_column)
    log(f"[mlst_f] Completed MLST summary for {report}")
