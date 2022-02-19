import os
import pandas as pd

# from snakemake.utils import validate

fastq_files = (
    pd.read_csv(config["fastq_files"], sep="\t", dtype={"sample_id": str, "fastq1": str, "fastq2": str})
    .set_index(["sample_id"], drop=False)
    .sort_index()
)

# validate(fastq_files, schema="../schemas/fastq_files.schema.yaml")

def get_final_output():
    # no trimming, use raw reads
    u = fastq_files[["fastq1", "fastq2"]].dropna()
    fastqs = u.fastq1.tolist() + u.fastq2.tolist()
    fastqs = [os.path.basename(fastq) for fastq in fastqs]
    fastqs = [fastq.rsplit('.fastq.gz', 1)[0] for fastq in fastqs]

    final_output = expand(
        "results/qc/fastqc/{fastq_file}_fastqc.html",
        fastq_file=fastqs,
    )
    return final_output

def get_fastqc_input(wildcards):
    u = fastq_files.dropna()
    u["index"] = range(u.shape[0])
    u = pd.wide_to_long(u, ["fastq"], i="index", j="read")
    u["basename_prefix"] = [os.path.basename(x).rstrip("fastq.gz") for x in u.fastq]
    u = u.set_index(["basename_prefix"], drop=False)
    input_file = u.loc[wildcards.fastq_prefix, ["fastq"]].tolist()
    return input_file
    