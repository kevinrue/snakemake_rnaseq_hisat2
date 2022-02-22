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
        else:
            shell("mv resources/genome.fa.gz resources/genome.fa")


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


## annotations


rule download_gene_annotations:
    input:
        ftp.remote(config["genes_gtf"])
    output:
        "resources/genes.gtf",
    log:
        "resources/genes.log",
    run:
        shell("mv {input} resources/genes.gtf.gz")
        if (bool(re.search(r".gz$", str(input)))):
            shell("gzip -d resources/genes.gtf.gz > {log} 2>&1")
        else:
            shell("mv resources/genes.gtf.gz resources/genes.gtf")
