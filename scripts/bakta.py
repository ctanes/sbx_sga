import sys
import os

def parse_file(filelines):
    if len(filelines) != 0:
        parsed_dict = {}
        for line in filelines:
            line = line.rstrip().split(":")
            try:
                key = line[0]
                value = line[1].strip()
            except:
                continue
            parsed_dict[key] = value
        return parsed_dict
    else:
        return {
            "CRISPR arrays": "NA",
            "hypotheticals": "NA",
            "tRNAs": "NA",
            "tmRNAs": "NA",
            "rRNAs": "NA",
            "CDSs": "NA",
            "N50": "NA",
            "Length": "NA",
            "GC": "NA",
        }


def get_annotation_stats(parsed_dict):
    crispr_count = parsed_dict["CRISPR arrays"]
    hypothetical_count = parsed_dict["hypotheticals"]
    trna_count = parsed_dict["tRNAs"]
    tmrna_count = parsed_dict["tmRNAs"]
    rrna_count = parsed_dict["rRNAs"]
    cds_count = parsed_dict["CDSs"]
    N50 = parsed_dict["N50"]
    genome_size = parsed_dict["Length"]
    gc = parsed_dict["GC"]
    return (
        crispr_count,
        hypothetical_count,
        trna_count,
        tmrna_count,
        rrna_count,
        cds_count,
        N50,
        genome_size,
        gc,
    )


def write_to_report(
    output,
    genome,
    crispr_count,
    hypothetical_count,
    trna_count,
    tmrna_count,
    rrna_count,
    cds_count,
    N50,
    genome_size,
    gc,
):
    sample = os.path.splitext(os.path.basename(genome))[0]
    with open(output, "w") as op:
        op.write(
            f"{sample}\t{genome_size}\t{cds_count}\t{N50}\t{rrna_count}\t{trna_count}\t{tmrna_count}\t{crispr_count}\t{hypothetical_count}\t{gc}\n"
        )
    # The file is closed after this function, and only the path is returned.
    return output


report = snakemake.input[0]
output = snakemake.output[0]

filelines = open(report, "r").readlines()
bakta_dict = parse_file(filelines)
(
    crispr_count,
    hypothetical_count,
    trna_count,
    tmrna_count,
    rrna_count,
    cds_count,
    N50,
    genome_size,
    gc,
) = get_annotation_stats(bakta_dict)
write_to_report(
    output,
    report,
    crispr_count,
    hypothetical_count,
    trna_count,
    tmrna_count,
    rrna_count,
    cds_count,
    N50,
    genome_size,
    gc,
)


