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


def test_parse():
    lines = ['Annotation:\n', 'test: 9\n', 'a: 10\n', '\n', 'd: 1343\n']
    assert parse_file(lines) == {"Annotation":"", "test":"9", "a":"10", "d":"1343"}

def write_to_report(output, genome, parsed_dict):
    sample = os.path.splitext(os.path.basename(genome))[0]
    with open(output, "w") as op:
        op.write(
            f"{sample}\t{parsed_dict['Length']}\t{parsed_dict['CDSs']}\t{parsed_dict['N50']}\t{parsed_dict['rRNAs']}\t{parsed_dict['tRNAs']}\t{parsed_dict['tmRNAs']}\t{parsed_dict['CRISPR arrays']}\t{parsed_dict['hypotheticals']}\t{parsed_dict['GC']}\n"
        )
    # The file is closed after this function, and only the path is returned.
    return output


