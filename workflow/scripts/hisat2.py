import sys
import os

# logging
sys.stderr = open(str(snakemake.log), "w")

fastq1 = ",".join(snakemake.input.fastq1)
fastq2 = ",".join(snakemake.input.fastq2)

os.popen(
    f"hisat2 \
        --threads {snakemake.params.threads} \
        -x {snakemake.input.index}/genome \
        -1 {fastq1} \
        -2 {fastq2} \
        {snakemake.params.hisat2} \
        --summary-file \
        {snakemake.log} \
    | samtools sort \
        -@ {snakemake.params.threads} \
        -o {snakemake.output} \
        - \
    && samtools index \
        {snakemake.output}"
)
