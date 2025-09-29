import os
import re
from typing import Callable, Iterable, Optional, Sequence, Set, Tuple


def _noop_log(message: str) -> None:
    """Default logger that ignores messages when no logger is provided."""
    pass


def open_report(
    report: str, log: Callable[[str], None] = _noop_log
) -> Tuple[str, Sequence[str]]:
    sample_name = os.path.basename(report.split("_sorted_winning.tab")[0])
    with open(report, "r") as report_obj:
        filelines = report_obj.readlines()
    top_lines = filelines[:20]
    log(
        f"[mash_f] Opened report {report} for sample {sample_name}; "
        f"total_lines={len(filelines)}, taking top {len(top_lines)}"
    )
    return sample_name, top_lines


def process_mash_line(line, log: Callable[[str], None] = _noop_log):
    try:
        parts = line.strip().split("\t")
        if len(parts) < 5:
            return None  # Malformed line

        identity = float(parts[0])
        median_multiplicity = float(parts[2])
        hits = int(parts[1].split("/")[0])

        # The last column(s) contain species info
        # Join remaining parts in case species field has tabs (some formats do)
        remainder = "\t".join(parts[4:])

        # Try matching common species name patterns
        species = extract_species_name(remainder)

        return species, median_multiplicity, identity, hits
    except Exception as e:
        log(f"[mash_f] Error processing line: {line}. Error: {e}")
        return "Unknown", 0, 0.0, 0


def extract_species_name(s):
    """Extract genus and species name from string, removing strain information."""
    # Format 1: ...-Staphylococcus_aureus.fna
    match = re.search(r"[-_](?P<species>[A-Z][a-z]+_[a-z]+(?:_[\w\d]+)?)\.fna", s)
    if match:
        full_name = match.group("species").replace("_", " ")
    else:
        # Format 2: ... Enterobacteria phage lambda, complete genome
        match = re.search(r"([A-Z][a-z]+(?: [a-z]+){1,3})", s)
        if match:
            full_name = match.group(1)
        else:
            return "Unknown"

    # Extract only genus and species (first two words)
    name_parts = full_name.split()
    if len(name_parts) >= 2:
        return f"{name_parts[0]} {name_parts[1]}"
    return full_name  # Return full name if it has fewer than 2 parts


def get_first_non_phage_hit(
    lines: Iterable[str], log: Callable[[str], None] = _noop_log
) -> Tuple[Optional[Tuple[str, float, float, int]], Optional[int]]:
    for idx, line in enumerate(lines):
        if "phage" not in line.lower():
            log(f"[mash_f] Found first non-phage hit at index {idx}")
            return process_mash_line(line, log), idx
    log("[mash_f] No non-phage hits detected in top lines")
    return None, None


def parse_report(
    top_lines: Sequence[str], log: Callable[[str], None] = _noop_log
) -> Set[str]:
    target_species: Set[str] = set()

    result = get_first_non_phage_hit(top_lines, log)

    if result == (None, None):
        log("[mash_f] parse_report returning empty set due to lack of non-phage hits")
        return target_species

    # Get top non-phage hit and its index
    (top_species, top_median_multiplicity, top_identity, top_hits), top_index = result

    if (top_identity >= 0.85) and (top_hits >= 100):
        target_species.add(top_species)
        log(f"[mash_f] Top hit passes thresholds: {top_species}")

    # Set the threshold for median multiplicity
    threshold = 0.05 * top_median_multiplicity
    log(f"[mash_f] Median multiplicity threshold set to {threshold}")

    # Iterate through the rest of the hits, excluding top_index
    for i, line in enumerate(top_lines):
        if i == top_index:
            continue
        species, median_multiplicity, identity, hits = process_mash_line(line, log)
        if (identity >= 0.85) and (hits >= 100):
            if any(term in species for term in ["phage", "Phage", "sp."]):
                continue
            if median_multiplicity >= threshold:
                target_species.add(species)
                log(f"[mash_f] Adding additional species {species}")

    log(f"[mash_f] parse_report returning {target_species}")
    return target_species


def contamination_call(target_set: Set[str], log: Callable[[str], None] = _noop_log):
    mash_dict = {}
    if len(target_set) <= 1:
        mash_dict["NA"] = ""
    else:
        species = " ".join(sorted(target_set))
        mash_dict["Contaminated"] = species
    log(f"[mash_f] contamination_call produced {mash_dict}")
    return mash_dict


def write_report(
    output: str,
    sample_name: str,
    mash_dict,
    log: Callable[[str], None] = _noop_log,
):
    # Expecting that the dictionary is just one key-value pair, so need to check that
    if len(mash_dict) == 1:
        status = list(mash_dict.keys())[0]
    else:
        # Raise error if dictionary is not proper length
        raise ValueError(
            f"Expected mash_dict to have exactly one key-value pair, but got {len(mash_dict)}."
        )
    log(f"[mash_f] Writing Mash report for {sample_name} with status={status}")
    with open(output, "w") as out:
        if status == "Contaminated":
            contaminated_spp = mash_dict[status]
            out.write(f"{sample_name}\tContaminated\t{contaminated_spp}\n")
        else:
            out.write(f"{sample_name}\tNA\tNA\n")
    log(f"[mash_f] Completed writing Mash report to {output}")
    return output
