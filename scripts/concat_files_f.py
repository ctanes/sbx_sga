import csv
import os
from typing import Iterable, TextIO


def log(message: str, log_file: TextIO) -> None:
    log_file.write(f"[concat_files_f] {message}\n")
    if hasattr(log_file, "flush"):
        log_file.flush()


def write_report(writer, report_reader: Iterable[list], sample_name: str, log_file: TextIO):
    log(f"Writing rows for sample {sample_name}", log_file)
    for row in report_reader:
        row.insert(0, sample_name)
        writer.writerow(row)


def summarize_all(
    report_paths,
    out_fp,
    folder_suffix="",
    header=True,
    log_file: TextIO = None,
):
    if log_file is None:
        raise ValueError("log_file must be provided for summarize_all")

    log(
        "Preparing to summarize "
        f"{len(report_paths)} reports to {out_fp} (header={header})",
        log_file,
    )
    first_non_empty = True
    header_first = []
    with open(out_fp, "w") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        for report_path in report_paths:
            log(f"Processing report {report_path}", log_file)
            if os.path.getsize(report_path) < 5:
                log(f"Skipping {report_path} because file is nearly empty", log_file)
                continue
            sample_name = os.path.basename(os.path.dirname(report_path))
            sample_name = sample_name.removesuffix(folder_suffix)
            with open(report_path, "r") as in_f:
                report_reader = csv.reader(in_f, delimiter="\t")
                if header:
                    header_line = next(report_reader)
                    header_line.insert(0, "SampleID")
                    if first_non_empty:
                        log(f"Setting header to {header_line}", log_file)
                        header_first = header_line
                        writer.writerow(header_first)
                    else:
                        if header_line != header_first:
                            log(
                                f"Header mismatch for sample {sample_name}: {header_line}",
                                log_file,
                            )
                            log(f"Expected header: {header_first}", log_file)
                            raise ValueError(
                                f"Headers in sample {sample_name} doesn't match the first file"
                            )
                write_report(writer, report_reader, sample_name, log_file)
            first_non_empty = False
    log(f"Finished writing combined report to {out_fp}", log_file)
