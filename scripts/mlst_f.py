import csv


def smush_column(line):
    line_parsed = []
    print(line)
    if line:
        schema = line[1]
        st = line[2]
        alleles = line[3:]
        alleles_joined = " ".join(alleles)
        line_parsed = [schema, st, alleles_joined]
    return line_parsed

def test_smush_column():
    mlst_row = ["marc.bacteremia.1020.a.fa", "saureus", "30", "arcC(2)", "aroE(2)", "glpF(2)", "yqiL(2)"]
    assert smush_column(mlst_row) == ["saureus", "30", "arcC(2) aroE(2) glpF(2) yqiL(2)"]


def write_to_report(report, output):
    with open(report, 'r') as f_in:
        reader = csv.reader(f_in, delimiter='\t')
        with open(output, "w") as op:
            writer = csv.writer(op, delimiter='\t')
            writer.writerow(["Schema", "ST", "Alleles"])
            for row in reader:
                smushed_column = smush_column(row)
                writer.writerow(smushed_column)
