import csv
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Callable


class ToolOutput(ABC):
    """Base class for tool outputs that enforces composability.
    Represents the output of a single tool for a single sample."""

    KEYS: list[str] = []

    def __init__(self, d: dict[str, str], keys: list[str]):
        self.name = "Abstract"
        assert all(key in d for key in keys), f"Missing keys: {keys} vs {d}"
        self.d = d

    @classmethod
    @abstractmethod
    def from_report(cls, fp: Path, log: Callable[[str], Any] = print) -> "ToolOutput":
        """Parses a report file and returns an instance of the tool output."""
        log(f"from_report not implemented for {cls.__name__}")
        return cls({}, [])


class BaktaOutput(ToolOutput):
    KEYS = [
        "Length",
        "GC",
        "N50",
        "CDSs",
        "tRNAs",
        "tmRNAs",
        "rRNAs",
        "hypotheticals",
        "CRISPR arrays",
    ]

    def __init__(self, d: dict[str, str]):
        super().__init__(d, self.KEYS)
        self.name = "Bakta"

    @classmethod
    def from_report(cls, fp: Path, log: Callable[[str], Any] = print) -> "BaktaOutput":
        log(f"Parsing Bakta report from {fp}")
        d = {}
        with open(fp) as f:
            for l in f:
                if len(l.split(":")) != 2:
                    continue
                key, value = l.split(":")
                key = key.strip()
                value = value.strip()

                if key in cls.KEYS:
                    d[key] = value

        return cls(d)


class CheckMOutput(ToolOutput):
    KEYS = ["Completeness", "Contamination"]

    def __init__(self, d: dict[str, str]):
        super().__init__(d, self.KEYS)
        self.name = "CheckM"

    @classmethod
    def from_report(cls, fp: Path, log: Callable[[str], Any] = print) -> "CheckMOutput":
        log(f"Parsing CheckM report from {fp}")
        d = {}
        with open(fp) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                for key in cls.KEYS:
                    if key in row:
                        d[key] = row[key]
                break  # Only read the first line
        return cls(d)


class MashOutput(ToolOutput):
    KEYS = ["Mash_Contamination", "Contaminated_Spp"]

    def __init__(self, d: dict[str, str]):
        super().__init__(d, self.KEYS)
        self.name = "Mash"

    @classmethod
    def from_report(cls, fp: Path, log: Callable[[str], Any] = print) -> "MashOutput":
        log(f"Parsing Mash report from {fp}")
        d = {}
        with open(fp) as f:
            pass
        return cls(d)


class MLSTOutput(ToolOutput):
    KEYS = ["Schema", "ST", "Alleles"]

    def __init__(self, d: dict[str, str]):
        super().__init__(d, self.KEYS)
        self.name = "MLST"

    @classmethod
    def from_report(cls, fp: Path, log: Callable[[str], Any] = print) -> "MLSTOutput":
        log(f"Parsing MLST report from {fp}")
        d = {}
        with open(fp) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample_name = row[0]
                d[cls.KEYS[0]] = row[1]
                d[cls.KEYS[1]] = row[2]
                d[cls.KEYS[2]] = " ".join(row[3:])
                break  # Only read the first line
        return cls(d)


class ShovillOutput(ToolOutput):
    KEYS = ["number of contigs", "min coverage", "max coverage", "mean coverage"]

    def __init__(self, d: dict[str, str]):
        super().__init__(d, self.KEYS)
        self.name = "Shovill"

    @staticmethod
    def _parse_header(h: str) -> dict[str, str]:
        hl = h.split()[1:]  # gets rid of the contig name
        header_dict = dict(item.split("=") for item in hl)
        return header_dict

    @classmethod
    def from_report(
        cls, fp: Path, log: Callable[[str], Any] = print
    ) -> "ShovillOutput":
        log(f"Parsing Shovill report from {fp}")
        headers = []
        with open(fp) as f:
            for l in f:
                if l.startswith(">"):
                    headers.append(l.strip())

        hs = [cls._parse_header(h) for h in headers]
        total_contigs = len(hs)
        min_cov = min(float(h["cov"]) for h in hs)
        max_cov = max(float(h["cov"]) for h in hs)
        total_length = sum(int(h["len"]) for h in hs)
        avg_cov = sum(int(h["len"]) * float(h["cov"]) for h in hs) / total_length
        rounded_cov = round(avg_cov, 2)

        return cls(
            {
                cls.KEYS[0]: str(total_contigs),
                cls.KEYS[1]: str(min_cov),
                cls.KEYS[2]: str(max_cov),
                cls.KEYS[3]: str(rounded_cov),
            }
        )


class SylphOutput(ToolOutput):
    KEYS = ["Taxonomic_abundance", "Contig_name"]

    def __init__(self, d: dict[str, str]):
        super().__init__(d, self.KEYS)
        self.name = "Sylph"

    @classmethod
    def from_report(cls, fp: Path, log: Callable[[str], Any] = print) -> "SylphOutput":
        log(f"Parsing Sylph report from {fp}")
        d = {}
        with open(fp) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                for key in cls.KEYS:
                    if key in row:
                        d[key] = row[key]
                break  # Only read the first line
        return cls(d)
