import os
import sys
from io import StringIO
import csv


def get_fasta_headers(f_in):
    headers = []
    for line in f_in:
        if line.startswith(">"):
            line = line.strip()
            headers.append(line)
    return headers


def test_get_fasta_headers():
    test_string = ">s1\nACG\nGCT\n>s2\nAC"
    test_stream = StringIO(test_string)
    headers = get_fasta_headers(test_stream)
    assert headers == [">s1", ">s2"]


def parse_header(header):
    header = header.split()[1:]  # gets rid of the contig name
    header_dict = dict(item.split("=") for item in header)
    header_dict["len"] = int(header_dict["len"])
    header_dict["cov"] = float(header_dict["cov"])
    keys_of_interest = ["len", "cov"]
    return {key: header_dict[key] for key in keys_of_interest if key in header_dict}


def test_parse_header():
    test_header = ">contig00001 len=122482 cov=39.3 corr=0 origname=Contig_88_39.3037_pilon sw=shovill-skesa/1.1.0 date=20250819"
    cov = parse_header(test_header)
    assert cov == {"len": 122482, "cov": 39.3}


def test_parse_header2():
    test_header = ">contig00001 len=122482 cov=39 corr=0 origname=Contig_88_39.3037_pilon sw=shovill-skesa/1.1.0 date=20250819"
    cov = parse_header(test_header)
    assert cov == {"len": 122482, "cov": 39}


def calc_cov_stats(contig_stats):
    total_ctgs = len(contig_stats)
    min_cov = min(c["cov"] for c in contig_stats)
    max_cov = max(c["cov"] for c in contig_stats)
    total_length = sum(c["len"] for c in contig_stats)
    avg_cov = sum(c["len"] * c["cov"] for c in contig_stats) / total_length

    rounded_cov = round(avg_cov, 2)
    return [total_ctgs, min_cov, max_cov, rounded_cov]


def test_calc_cov_stats():
    d = [{"len": 2, "cov": 10}, {"len": 5, "cov": 2}]
    assert calc_cov_stats(d) == [2, 2, 10, 4.29]


def write_shovill_stats(fp_in, fp_out):
    with open(fp_in, "r") as f_in:
        headers = get_fasta_headers(f_in)

    with open(fp_out, "w") as op:
        writer = csv.writer(op, delimiter="\t")
        writer.writerow(
            ["number of contigs", "min coverage", "max coverage", "mean coverage"]
        )
        if headers:
            contig_stats = []
            for header in headers:
                contig_stats.append(parse_header(header))
            final_stats = calc_cov_stats(contig_stats)
            writer.writerow(final_stats)
        else:
            writer.writerow([0, 0, 0, 0])
