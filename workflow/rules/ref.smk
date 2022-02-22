## genome


rule download_genome:
    input:
        ftp.remote(config["genome_fasta"])
    output:
        "resources/genome.fa",
    log:
        "resources/genome.log",
    run:
        shell("mv {input} resources/genome.fa.gz")
        if (bool(re.search(r".gz$", str(input)))):
            shell("gzip -d resources/genome.fa.gz > {log} 2>&1")


## HISAT2

rule index_genome:
    input:
        "resources/genome.fa"
    output:
        directory("resources/index_genome")
    log:
        "resources/index_genome.log"
    conda:
        "../envs/hisat2.yaml"
    shell:
        """
        mkdir -p {output}
        hisat2-build {input} {output}/genome {log} 2>&1
        """
