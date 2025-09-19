ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
try:
    SBX_SGA_VERSION = get_ext_version("sbx_sga")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_SGA_VERSION = "0.0.0"


localrules:
    all_sga_virus,


rule all_sga_virus:
    input:
        expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_virus_summary.tsv",
            sample=Samples,
        ),
        f"{ISOLATE_FP}/reports/genomad_plasmid_classification.report",
        f"{ISOLATE_FP}/reports/genomad_virus_classification.report",
        f"{ISOLATE_FP}/reports/genomad_plasmid_genes.report",
        f"{ISOLATE_FP}/reports/genomad_virus_genes.report",


rule sga_genomad_download_db:
    """Download Genomad database"""
    output:
        version=Path(Cfg["sbx_sga"]["genomad_ref"]) / "genomad_db" / "version.txt",
    log:
        LOG_FP / "sga_genomad_download_db.log",
    benchmark:
        BENCHMARK_FP / "sga_genomad_download_db.tsv"
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        GENOMAD_DB_DIR=$(dirname {output.version})
        genomad download-database $(dirname "$GENOMAD_DB_DIR") > {log} 2>&1 || true
        """


rule sga_genomad_end_to_end:
    """Run Genomad end-to-end pipeline"""
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
        db_version=Path(Cfg["sbx_sga"]["genomad_ref"]) / "genomad_db" / "version.txt",
    output:
        plasmid_summary=ISOLATE_FP
        / "genomad"
        / "{sample}"
        / "{sample}_summary"
        / "{sample}_plasmid_summary.tsv",
        virus_summary=ISOLATE_FP
        / "genomad"
        / "{sample}"
        / "{sample}_summary"
        / "{sample}_virus_summary.tsv",
        plasmid_gene_summary=ISOLATE_FP
        / "genomad"
        / "{sample}"
        / "{sample}_summary"
        / "{sample}_plasmid_genes.tsv",
        virus_gene_summary=ISOLATE_FP
        / "genomad"
        / "{sample}"
        / "{sample}_summary"
        / "{sample}_virus_genes.tsv",
    log:
        LOG_FP / "genomad_end_to_end_{sample}.log",
    benchmark:
        BENCHMARK_FP / "genomad_end_to_end_{sample}.tsv"
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        ASSEMBLY_SUMMARY_DIR=$(dirname {output.plasmid_summary})
        DB_DIR=$(dirname {input.db_version})
        
        if [ ! -s {input.contigs} ]; then
            mkdir -p ${{ASSEMBLY_SUMMARY_DIR}}
            touch {output.plasmid_summary}
            touch {output.virus_summary}
            touch {output.plasmid_gene_summary}
            touch {output.virus_gene_summary}
        else
            genomad end-to-end --cleanup --splits 8 {input.contigs} $(dirname "${{ASSEMBLY_SUMMARY_DIR}}") $DB_DIR > {log} 2>&1
        fi
        """


rule sga_plasmid_classification_summary:
    """Summarizes the plasmid classifications per sample"""
    input:
        reports=expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_plasmid_summary.tsv",
            sample=Samples,
        ),
    output:
        plasmid_summary_report=ISOLATE_FP
        / "reports"
        / "genomad_plasmid_classification.report",
    params:
        suffix="_summary",
        header=True,
    log:
        LOG_FP / "sga_plasmid_classification_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_plasmid_classification_summary.tsv"
    script:
        "scripts/concat_files.py"


rule sga_phage_classification_summary:
    """Summarizes the phage classifications per sample"""
    input:
        reports=expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_virus_summary.tsv",
            sample=Samples,
        ),
    output:
        virus_summary_report=ISOLATE_FP
        / "reports"
        / "genomad_virus_classification.report",
    params:
        suffix="_summary",
        header=True,
    log:
        LOG_FP / "sga_phage_classification_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_phage_classification_summary.tsv"
    script:
        "scripts/concat_files.py"


rule sga_plasmid_gene_summary:
    """Summarizes the genes associated with each of the plasmid classifications per sample"""
    input:
        reports=expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_virus_genes.tsv",
            sample=Samples,
        ),
    output:
        plasmid_gene_report=ISOLATE_FP / "reports" / "genomad_plasmid_genes.report",
    params:
        suffix="_summary",
        header=True,
    log:
        LOG_FP / "sga_plasmid_gene_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_plasmid_gene_summary.tsv"
    script:
        "scripts/concat_files.py"


rule sga_virus_gene_summary:
    """Summarizes the genes associated with each of the phage classifications per sample"""
    input:
        reports=expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_plasmid_genes.tsv",
            sample=Samples,
        ),
    output:
        virus_gene_report=ISOLATE_FP / "reports" / "genomad_virus_genes.report",
    params:
        suffix="_summary",
        header=True,
    log:
        LOG_FP / "sga_virus_gene_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_virus_gene_summary.tsv"
    script:
        "scripts/concat_files.py"
