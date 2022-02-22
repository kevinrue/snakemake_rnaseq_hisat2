## samtools


rule idxstats_on_bam:
    input:
        "results/hisat2/{sample_id}.bam",
    output:
        "results/qc/samtools/idxstats/{sample_id}",
    log:
        "results/qc/samtools/idxstats/{sample_id}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats \
            {input} \
            > {output} \
            2> {log}
        """


rule flagstat_on_bam:
    input:
        "results/hisat2/{sample_id}.bam",
    output:
        "results/qc/samtools/flagstat/{sample_id}",
    log:
        "results/qc/samtools/flagstat/{sample_id}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat \
            {input} \
            > {output} \
            2> {log}
        """


rule picard_alignment_metrics_on_bam:
    input:
        bam="results/hisat2/{sample_id}.bam",
        genome="resources/genome.fa",
    output:
        "results/qc/picard/CollectAlignmentSummaryMetrics/{sample_id}",
    log:
        "results/qc/picard/CollectAlignmentSummaryMetrics/{sample_id}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard CollectAlignmentSummaryMetrics \
            -I {input.bam} \
            -O {output} \
            -R {input.genome} \
            2> {log} \
            2>&1
        """


rule picard_insert_size_on_bam:
    input:
        bam="results/hisat2/{sample_id}.bam",
        genome="resources/genome.fa",
    output:
        "results/qc/picard/CollectInsertSizeMetrics/{sample_id}",
    log:
        "results/qc/picard/CollectInsertSizeMetrics/{sample_id}.log",
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard CollectInsertSizeMetrics \
            -I {input.bam} \
            -O {output} \
            -R {input.genome} \
            -H {output}.pdf \
            > {log} \
            2>&1
        """
