## FastQC


rule fastqc_on_fastq:
    input:
        get_fastqc_input,
    output:
        fastqc="results/qc/fastqc/{fastq_prefix}_fastqc.html",
    log:
        "results/qc/fastqc/{fastq_prefix}_fastqc.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc \
            -o results/qc/fastqc \
            --nogroup \
            {input} \
            > {log} \
            2>&1
        """


## MultiQC


rule multiqc_on_reads:
    input:
        get_reads_reports(),
    output:
        "results/reports/multiqc/reads.html",
    log:
        "results/reports/multiqc/reads.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc \
            -n reads.html \
            -o results/reports/multiqc \
            results/qc/fastqc \
            > {log} \
            2>&1
        """
