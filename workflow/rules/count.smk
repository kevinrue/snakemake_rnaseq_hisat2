## featurecounts


rule featurecounts_on_bam:
    input:
        bam=get_featurecounts_input(),
        gtf="resources/genes.gtf.gz",
    output:
        "results/featurecounts/counts",
    params:
        threads=config["featurecounts"]["threads"],
        featurecounts=config["featurecounts"]["options"],
    log:
        "results/featurecounts/counts.log",
    conda:
        "../envs/subread.yaml"
    shell:
        """
        featureCounts \
            -a {input.gtf} \
            -o {output} \
            -T {params.threads} \
            {params.featurecounts} \
           {input.bam} \
            > {log} \
            2>&1
        """
