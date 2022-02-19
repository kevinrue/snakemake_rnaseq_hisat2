

rule fastqc_on_fastq:
    input:
        get_fastqc_input,
    output:
        fastqc="results/qc/fastqc/{fastq_prefix}_fastqc.html",
    log:
        "results/qc/fastqc/{fastq_prefix}_fastqc.html.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc \
            -o results/qc/fastqc \
            --nogroup \
            {input} \
            > results/qc/fastqc/{wildcards.fastq_prefix}_fastqc.html.log \
            2>&1
        """