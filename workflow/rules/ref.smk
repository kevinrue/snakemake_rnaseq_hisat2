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
