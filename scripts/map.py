import pandas as pd


assembly_qc = {
    "shovill": {
        "Total_contigs": "contig_count",
        "Min_coverage": "min_contig_coverage",
        "Max_coverage": "max_contig_coverage",
        "Total_length": "genome_size",
        "Average_coverage": "avg_contig_coverage",
    },
    "bakta": {
        "GC": "gc_content",
        "N50": "n50",
        "CDSs": "cds",
    },
    "checkm": {
        "Completeness": "completeness",
        "Contamination": "contamination",
    },
}

taxonomic_assignment = {
    "mlst": {
        "classification": "classification",
        "allele_assignment": "comment",
    },
    "sylph": {
        "Contig_name": "classification",
    },
}

contaminant = {
    "mash": {
        "hits_per_thousand": "confidence",
        "species": "classification",
    }
}

antimicrobial = {
    "abritamr": {
        "Contig id": "contig_id",
        "Gene symbol": "gene_symbol",
        "Sequence name": "gene_name",
        "Accession of closest sequence": "accession",
        "Element type": "element_type",
        "Subclass": "resistance_product",
    }
}

virus = {
    "genomad_plasmid_summary": {
        "seq_name": "contig_id",
    },
    "genomad_virus_summary": {
        # TODO: Sometimes this field contains extra stuff after a "|" character
        "seq_name": "contig_id",
    },
    "genomad_plasmid_genes": {
        # TODO: Figure out how to include transforms here maybe?? (split this field on "_")
        "contig_id": "gene",
        "gene_name": "gene",
    },
    "genomad_virus_genes": {
        "contig_id": "gene",
        "gene_name": "gene",
    },
}


def reduce_dataframe(df: pd.DataFrame, tool: str) -> pd.DataFrame:
    """Reduce dataframe to only relevant columns based on tool type and rename for consistency with database model."""
    if tool in assembly_qc:
        relevant_keys = assembly_qc[tool]
    elif tool in taxonomic_assignment:
        relevant_keys = taxonomic_assignment[tool]
    elif tool in contaminant:
        relevant_keys = contaminant[tool]
    elif tool in antimicrobial:
        relevant_keys = antimicrobial[tool]
    elif tool in virus:
        relevant_keys = virus[tool]
    else:
        raise ValueError(f"Unknown tool: {tool}")

    relevant_columns = [
        col for col in df.columns if col in relevant_keys or col == "SampleID"
    ]
    df = df[relevant_columns]
    df = df.rename(columns={k: v for k, v in relevant_keys.items()})
    return df
