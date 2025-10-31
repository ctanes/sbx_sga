ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
TOOLS = {
    "shovill": ["number of contigs", "min coverage", "max coverage", "mean coverage"],
    "checkm": ["Completeness", "Contamination"],
    "sylph": ["Taxonomic_abundance", "Contig_name"],
    "mlst": ["Schema", "ST", "Alleles"],
    "bakta": [
        "Length",
        "GC",
        "N50",
        "CDSs",
        "tRNAs",
        "tmRNAs",
        "rRNAs",
        "hypotheticals",
        "CRISPR arrays",
    ],
    "mash": ["Mash_Contamination", "Contaminated_Spp"],
}

try:
    SBX_SGA_VERSION = get_ext_version("sbx_sga")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_SGA_VERSION = "0.0.0"


localrules:
    all_sga,
    sga_report,


rule all_sga:
    input:
        expand(ISOLATE_FP / "quast" / "{sample}" / "report.tsv", sample=Samples),
        f"{ISOLATE_FP}/final_summary.tsv",
        f"{ISOLATE_FP}/reports/amr.report",


## Assembly
rule sga_shovill:
    input:
        rp1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        rp2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    log:
        LOG_FP / "sga_shovill_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_shovill_{sample}.tsv"
    conda:
        "envs/shovill.yml"
    resources:
        mem_mb=32000,
    shell:
        """
        (shovill --force --assembler skesa --outdir $(dirname {output.contigs}) --R1 {input.rp1}  --R2 {input.rp2} &> {log} &&
        mv $(dirname {output.contigs})/contigs.fa {output.contigs}) ||
        touch {output.contigs}
        """


# Taxonomic classification
rule sga_sylph:
    input:
        rp1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        rp2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        report=ISOLATE_FP / "sylph" / "{sample}" / "{sample}.tsv",
    threads: 8
    params:
        ref=Cfg["sbx_sga"]["sylph_ref"],
    log:
        LOG_FP / "sga_sylph_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_sylph_{sample}.tsv"
    conda:
        "envs/sylph.yml"
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        """
    if [ $(zcat {input.rp1} | wc -l) -ge 4 ] && [ $(zcat {input.rp2} | wc -l) -ge 4 ]; then
        sylph sketch -1 "{input.rp1}" -2 "{input.rp2}" -t 1 -d "$(dirname "{output.report}")" > "{log}" 2>&1
        sylph profile "{params.ref}" "$(dirname "{output.report}")"/*.sylsp -t {threads} -o "{output.report}" >> "{log}" 2>&1
    else
        touch {output.report}
    fi
        """


### Assembly QC
rule sga_checkm:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        quality_report=ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv",
    params:
        ref=Cfg["sbx_sga"]["checkm_ref"],
    log:
        LOG_FP / "sga_checkm_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_checkm_{sample}.tsv"
    conda:
        "envs/checkm2.yml"
    resources:
        mem_mb=16000,
    shell:
        """
        if [ -s {input.contigs} ]; then
            checkm2 predict \\
            --force \\
            -x fa \\
            -i {input.contigs} \\
            -o $(dirname {output.quality_report}) \\
            --database_path {params.ref} \\
            &> {log} || touch {output.quality_report}
        else
            touch {output.quality_report}
        fi
        """


rule sga_quast:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        quast_dir=ISOLATE_FP / "quast" / "{sample}" / "report.tsv",
    log:
        LOG_FP / "sga_quast_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_quast_{sample}.tsv"
    conda:
        "envs/quast.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            quast.py \\
                -o $(dirname {output.quast_dir}) \\
                {input.contigs} \\
                &> {log} || touch {output.quast_dir}
        else
            touch {output.quast_dir}
        fi
        """


rule sga_mash:
    input:
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        agg=temp(ISOLATE_FP / "mash" / "{sample}" / "{sample}.fastq"),
        win=temp(ISOLATE_FP / "mash" / "{sample}" / "{sample}_winning.tab"),
        sort=ISOLATE_FP / "mash" / "{sample}" / "{sample}_sorted_winning.tab",
    params:
        ref=Cfg["sbx_sga"]["mash_ref"],
    log:
        LOG_FP / "sga_mash_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mash_{sample}.tsv"
    conda:
        "envs/mash.yml"
    resources:
        mem_mb=16000,
    shell:
        """
        zcat {input.reads} > {output.agg}

        if [ -s {output.agg} ]; then
            mash screen -w -p 8 {params.ref} {output.agg} > {output.win} 2> {log}
            sort -gr {output.win} > {output.sort} 2>> {log}
        else
            touch {output.win} {output.sort}
        fi
        """


# Typing
rule sga_mlst:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        mlst=ISOLATE_FP / "mlst" / "{sample}" / "{sample}.mlst",
    log:
        LOG_FP / "sga_mlst_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mlst_{sample}.tsv"
    conda:
        "envs/mlst.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            mlst --nopath {input.contigs} > {output.mlst} 2> {log}
        else
            mkdir -p $(dirname {output.mlst})
            touch {output.mlst}
        fi
        """


### Annotation
rule sga_bakta:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        bakta=ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt",
    params:
        ref=Cfg["sbx_sga"]["bakta_ref"],
    log:
        LOG_FP / "sga_bakta_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_bakta_{sample}.tsv"
    conda:
        "envs/bakta.yml"
    resources:
        mem_mb=32000,
        runtime=180,
    shell:
        """
        if [ -s {input.contigs} ]; then
            bakta --force --db {params.ref} \\
            --output $(dirname {output.bakta}) \\
            --prefix {wildcards.sample} \\
            --skip-plot {input.contigs} \\
            &> {log}
        else
            touch {output.bakta}
        fi
        """


### AMR Profiling
rule sga_abritamr:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        abritamr=ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out",
    log:
        LOG_FP / "sga_abritamr_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_abritamr_{sample}.tsv"
    conda:
        "envs/abritamr.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            abritamr run \\
                --contigs {input.contigs} \\
                --prefix {wildcards.sample} \\
                &> {log}
            mv {wildcards.sample} $(dirname $(dirname {output.abritamr}))
        else
            mkdir -p $(dirname {output.abritamr})
            touch {output.abritamr}
        fi     
    """


rule sga_report:
    input:
        shovill=expand(ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa", sample=Samples),
        sylph=expand(ISOLATE_FP / "sylph" / "{sample}" / "{sample}.tsv", sample=Samples),
        checkm=expand(ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv", sample=Samples),
        mlst=expand(ISOLATE_FP / "mlst" / "{sample}" / "parsed_mlst.txt", sample=Samples),
        bakta=expand(ISOLATE_FP / "bakta" / "{sample}" / "parsed_summary.txt", sample=Samples),
        mash=expand(ISOLATE_FP / "mash" / "{sample}" / "{sample}_summary.tsv", sample=Samples),
        abritamr=expand(ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out", sample=Samples),
    output:
        tool_reports=expand(ISOLATE_FP / "reports" / "{tool}.tsv", tool=[
            "shovill",
            "sylph",
            "checkm",
            "mlst",
            "bakta",
            "mash",
            "abritamr",
        ]),
        assembly_qcs=ISOLATE_FP / "assembly_qcs.tsv",
        taxonomic_assignments=ISOLATE_FP / "taxonomic_assignments.tsv",
        antimicrobials=ISOLATE_FP / "antimicrobials.tsv",
    log:
        LOG_FP / "sga_report.log",
    benchmark:
        BENCHMARK_FP / "sga_report.tsv"
    script:
        "scripts/summarize_all.py"
