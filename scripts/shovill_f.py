import csv
import os
from io import StringIO
from typing import Callable, Iterable, List


def _noop_log(message: str) -> None:
    """Default logger that ignores messages when no logger is provided."""
    pass


def get_fasta_headers(
    f_in: Iterable[str], log: Callable[[str], None] = _noop_log
) -> List[str]:
    headers = []
    for line in f_in:
        if line.startswith(">"):
            line = line.strip()
            headers.append(line)
            log(f"[shovill_f] Found header: {line}")
    log(f"[shovill_f] Total headers parsed: {len(headers)}")
    return headers


def test_get_fasta_headers():
    test_string = ">s1\nACG\nGCT\n>s2\nAC"
    test_stream = StringIO(test_string)
    buffer = StringIO()

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    headers = get_fasta_headers(test_stream, test_log)
    assert headers == [">s1", ">s2"]


def parse_header(header: str, log: Callable[[str], None] = _noop_log):
    header = header.split()[1:]  # gets rid of the contig name
    header_dict = dict(item.split("=") for item in header)
    header_dict["len"] = int(header_dict["len"])
    header_dict["cov"] = float(header_dict["cov"])
    keys_of_interest = ["len", "cov"]
    parsed = {key: header_dict[key] for key in keys_of_interest if key in header_dict}
    log(f"[shovill_f] Parsed header into {parsed}")
    return parsed


def test_parse_header():
    test_header = ">contig00001 len=122482 cov=39.3 corr=0 origname=Contig_88_39.3037_pilon sw=shovill-skesa/1.1.0 date=20250819"
    buffer = StringIO()

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    cov = parse_header(test_header, test_log)
    assert cov == {"len": 122482, "cov": 39.3}


def test_parse_header2():
    test_header = ">contig00001 len=122482 cov=39 corr=0 origname=Contig_88_39.3037_pilon sw=shovill-skesa/1.1.0 date=20250819"
    buffer = StringIO()

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    cov = parse_header(test_header, test_log)
    assert cov == {"len": 122482, "cov": 39}


def calc_cov_stats(contig_stats, log: Callable[[str], None] = _noop_log):
    total_ctgs = len(contig_stats)
    min_cov = min(c["cov"] for c in contig_stats)
    max_cov = max(c["cov"] for c in contig_stats)
    total_length = sum(c["len"] for c in contig_stats)
    avg_cov = sum(c["len"] * c["cov"] for c in contig_stats) / total_length

    rounded_cov = round(avg_cov, 2)
    stats = [total_ctgs, min_cov, max_cov, rounded_cov]
    log(
        "[shovill_f] Calculated coverage stats: "
        f"total_contigs={total_ctgs}, min={min_cov}, max={max_cov}, mean={rounded_cov}"
    )
    return stats


def test_calc_cov_stats():
    d = [{"len": 2, "cov": 10}, {"len": 5, "cov": 2}]
    buffer = StringIO()

    def test_log(message: str) -> None:
        buffer.write(f"{message}\n")

    assert calc_cov_stats(d, test_log) == [2, 2, 10, 4.29]


def write_shovill_stats(
    fp_in: str, fp_out: str, log: Callable[[str], None] = _noop_log
):
    log(f"[shovill_f] Writing Shovill statistics from {fp_in} to {fp_out}")
    if not os.path.exists(fp_in):
        raise FileNotFoundError(f"Input FASTA {fp_in} does not exist")
    with open(fp_in, "r") as f_in:
        headers = get_fasta_headers(f_in, log)

    with open(fp_out, "w") as op:
        writer = csv.writer(op, delimiter="\t")
        writer.writerow(
            ["number of contigs", "min coverage", "max coverage", "mean coverage"]
        )
        if headers:
            contig_stats = []
            for header in headers:
                contig_stats.append(parse_header(header, log))
            final_stats = calc_cov_stats(contig_stats, log)
            writer.writerow(final_stats)
            log(f"[shovill_f] Wrote stats row: {final_stats}")
        else:
            writer.writerow([0, 0, 0, 0])
            log("[shovill_f] No headers found; wrote zero stats")
