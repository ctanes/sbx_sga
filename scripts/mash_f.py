import os
import re
from typing import Iterable, Optional, Sequence, Set, TextIO, Tuple


def log(message: str, log_file: TextIO) -> None:
    log_file.write(f"[mash_f] {message}\n")
    if hasattr(log_file, "flush"):
        log_file.flush()


def open_report(report: str, log_file: TextIO) -> Tuple[str, Sequence[str]]:
    sample_name = os.path.basename(report.split("_sorted_winning.tab")[0])
    with open(report, "r") as report_obj:
        filelines = report_obj.readlines()
    top_lines = filelines[:20]
    log(
        f"Opened report {report} for sample {sample_name}; "
        f"total_lines={len(filelines)}, taking top {len(top_lines)}",
        log_file,
    )
    return sample_name, top_lines


def process_mash_line(line: str, log_file: TextIO) -> Tuple[str, float, float, int]:
    line_list = line.rstrip().split("\t")
    species_line = line_list[-1]
    matches = re.findall(r"N[A-Z]_[0-9A-Z]+\.[0-9]", species_line)
    if not matches:
        try:
            matches = re.findall(r"[A-Z]{2}_[0-9]+\.[0-9]", species_line)
        except Exception as exc:
            raise ValueError(f"No match found in species_line: {species_line}") from exc
    split_char = matches[0]
    species_split = line.split(split_char)[1].lstrip()
    species = " ".join(species_split.split()[:2])
    median_multiplicity = float(line_list[2])
    identity = float(line_list[0])
    hits = int(line_list[1].split("/")[0])
    log(
        "Processed mash line with species="
        f"{species}, identity={identity}, hits={hits}, median_multiplicity={median_multiplicity}",
        log_file,
    )
    return species, median_multiplicity, identity, hits


def get_first_non_phage_hit(lines: Iterable[str], log_file: TextIO) -> Tuple[Optional[Tuple[str, float, float, int]], Optional[int]]:
    for idx, line in enumerate(lines):
        if "phage" not in line.lower():
            log(f"Found first non-phage hit at index {idx}", log_file)
            return process_mash_line(line, log_file), idx
    log("No non-phage hits detected in top lines", log_file)
    return None, None


def parse_report(top_lines: Sequence[str], log_file: TextIO) -> Set[str]:
    target_species: Set[str] = set()

    result = get_first_non_phage_hit(top_lines, log_file)

    if result == (None, None):
        log("parse_report returning empty set due to lack of non-phage hits", log_file)
        return target_species

    # Get top non-phage hit and its index
    (top_species, top_median_multiplicity, top_identity, top_hits), top_index = result

    if (top_identity >= 0.85) and (top_hits >= 100):
        target_species.add(top_species)
        log(f"Top hit passes thresholds: {top_species}", log_file)

    # Set the threshold for median multiplicity
    threshold = 0.05 * top_median_multiplicity
    log(f"Median multiplicity threshold set to {threshold}", log_file)

    # Iterate through the rest of the hits, excluding top_index
    for i, line in enumerate(top_lines):
        if i == top_index:
            continue
        species, median_multiplicity, identity, hits = process_mash_line(line, log_file)
        if (identity >= 0.85) and (hits >= 100):
            if any(term in species for term in ["phage", "Phage", "sp."]):
                continue
            if median_multiplicity >= threshold:
                target_species.add(species)
                log(f"Adding additional species {species}", log_file)

    log(f"parse_report returning {target_species}", log_file)
    return target_species


def contamination_call(target_set: Set[str], log_file: TextIO):
    mash_dict = {}
    if len(target_set) <= 1:
        mash_dict["NA"] = ""
    else:
        species = " ".join(sorted(target_set))
        mash_dict["Contaminated"] = species
    log(f"contamination_call produced {mash_dict}", log_file)
    return mash_dict


def write_report(output: str, sample_name: str, mash_dict, log_file: TextIO):
    # Expecting that the dictionary is just one key-value pair, so need to check that
    if len(mash_dict) == 1:
        status = list(mash_dict.keys())[0]
    else:
        # Raise error if dictionary is not proper length
        raise ValueError(
            f"Expected mash_dict to have exactly one key-value pair, but got {len(mash_dict)}."
        )
    log(f"Writing Mash report for {sample_name} with status={status}", log_file)
    with open(output, "w") as out:
        if status == "Contaminated":
            contaminated_spp = mash_dict[status]
            out.write(f"{sample_name}\tContaminated\t{contaminated_spp}\n")
        else:
            out.write(f"{sample_name}\tNA\tNA\n")
    log(f"Completed writing Mash report to {output}", log_file)
    return output
