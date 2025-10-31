ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
try:
    SBX_SGA_VERSION = get_ext_version("sbx_sga")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_SGA_VERSION = "0.0.0"


virus_outputs = [
    "plasmid_summary",
    "virus_summary",
    "plasmid_genes",
    "virus_genes",
]


localrules:
    all_sga_virus,
    sga_virus_report,


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
        expand(
            ISOLATE_FP / "reports" / f"genomad_{out}.tsv",
            out=virus_outputs,
        ),


rule sga_genomad_end_to_end:
    """Run Genomad end-to-end pipeline"""
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
        db_version=Path(Cfg["sbx_sga"]["genomad_ref"]) / "genomad_db" / "version.txt",
    output:
        expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_{out}.tsv",
            out=virus_outputs,
        ),
    log:
        LOG_FP / "genomad_end_to_end_{sample}.log",
    benchmark:
        BENCHMARK_FP / "genomad_end_to_end_{sample}.tsv"
    conda:
        "envs/sga_virus.yml"
    shell:
        """
        ASSEMBLY_SUMMARY_DIR=$(dirname {output[0]})
        DB_DIR=$(dirname {input.db_version})
        
        if [ ! -s {input.contigs} ]; then
            echo "Empty input assembly file {input.contigs}, creating empty summary files" > {log}
            mkdir -p ${{ASSEMBLY_SUMMARY_DIR}}
            for out_file in {output}; do
                touch ${{out_file}}
            done
        else
            genomad end-to-end --cleanup --splits 8 {input.contigs} $(dirname "${{ASSEMBLY_SUMMARY_DIR}}") $DB_DIR > {log} 2>&1
        fi
        """


rule sga_virus_report:
    input:
        expand(
            ISOLATE_FP
            / "genomad"
            / "{sample}"
            / "{sample}_summary"
            / "{sample}_{out}.tsv",
            sample=Samples,
            out=virus_outputs,
        ),
    output:
        tool_reports=expand(
            ISOLATE_FP / "reports" / f"genomad_{out}.tsv", out=virus_outputs
        ),
        virus=ISOLATE_FP / "virus.tsv",
    log:
        LOG_FP / "sga_virus_report.log",
    benchmark:
        BENCHMARK_FP / "sga_virus_report.tsv"
    script:
        "scripts/summarize_virus.py"
