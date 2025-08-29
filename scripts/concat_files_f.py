import os
import csv

def write_report(writer, report_reader, sample_name):
    for row in report_reader:
        row.insert(0, sample_name)
        writer.writerow(row)

def summarize_all(report_paths, out_fp, folder_suffix="", header=True):
    first_non_empty = True
    header_first = []
    with open(out_fp, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        for report_path in report_paths:
            if os.path.getsize(report_path) < 5:
                continue
            sample_name = os.path.basename(os.path.dirname(report_path))
            sample_name = sample_name.removesuffix(folder_suffix)
            with open(report_path, "r") as in_f:
                report_reader = csv.reader(in_f, delimiter='\t')
                if header:
                    header_line = next(report_reader)
                    header_line.insert(0, "SampleID")
                    if first_non_empty:
                        header_first = header_line
                        writer.writerow(header_first)
                    else:
                        if (header_line != header_first):
                            print(header_line)
                            print(header_first)
                            raise ValueError(f"Headers in sample {sample_name} doesn't match the first file")
                write_report(writer, report_reader, sample_name)
            first_non_empty = False
