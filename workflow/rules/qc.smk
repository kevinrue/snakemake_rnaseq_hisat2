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


rule multiqc_on_fastqc:
    input:
        get_fastqc_reports(),
    output:
        "results/reports/multiqc/fastq.html",
    log:
        "results/reports/multiqc/fastq.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc \
            -n fastq.html \
            -o results/reports/multiqc \
            results/qc/fastqc \
            > {log} \
            2>&1
        """
