## HISAT2


rule hisat2_on_fastq:
    input:
        unpack(get_hisat2_input),
    output:
        "results/hisat2/{sample_id}.bam",
    params:
        hisat2=config["hisat2"]["options"],
        threads=config["hisat2"]["threads"],
    log:
        "results/hisat2/{sample_id}.log",
    conda:
        "../envs/hisat2.yaml"
    script:
        "../scripts/hisat2.py"
