from io import StringIO
from shovill_f import (
    get_fasta_headers,
    parse_header,
    calc_cov_stats,
    write_shovill_stats,
)


def test_get_fasta_headers():
    fasta = ">contig1 len=10 cov=20\nACGT\n>contig2 len=5 cov=30\nTTTT\n"
    headers = get_fasta_headers(StringIO(fasta))
    assert headers == [">contig1 len=10 cov=20", ">contig2 len=5 cov=30"]


def test_parse_header():
    header = (
        ">contig00001 len=122482 cov=39.3 corr=0 origname=Contig_88_39.3037_pilon "
        "sw=shovill-skesa/1.1.0 date=20250819"
    )
    assert parse_header(header) == {"len": 122482, "cov": 39.3}


def test_calc_cov_stats():
    contigs = [{"len": 2, "cov": 10}, {"len": 5, "cov": 2}]
    assert calc_cov_stats(contigs) == [2, 2, 10, 4.29]


def test_write_shovill_stats(tmp_path):
    fasta = (
        ">contig00001 len=129214 cov=78.8\nAAA\n"
        ">contig00039 len=24105 cov=86.0\nTTT\n"
    )
    fp_in = tmp_path / "input.fa"
    fp_in.write_text(fasta)
    fp_out = tmp_path / "out.tsv"
    write_shovill_stats(fp_in, fp_out)

    assert fp_out.exists()
    lines = fp_out.read_text().splitlines()
    assert lines[0] == "number of contigs\tmin coverage\tmax coverage\tmean coverage"
    assert lines[1] == "2\t78.8\t86.0\t79.93"
