import os
import re
import pandas as pd

from snakemake.remote import FTP

ftp = FTP.RemoteProvider()

# from snakemake.utils import validate

# validate(config, schema="../schemas/config.schema.yaml")

fastq_files = (
    pd.read_csv(
        config["fastq_files"],
        sep="\t",
        dtype={"sample_id": str, "fastq1": str, "fastq2": str},
    )
    .set_index(["sample_id"], drop=False)
    .sort_index()
)

# validate(fastq_files, schema="../schemas/fastq_files.schema.yaml")


def get_final_output():
    # no trimming, use raw reads
    u = fastq_files[["fastq1", "fastq2"]].dropna()
    fastqs = u.fastq1.tolist() + u.fastq2.tolist()
    fastqs = [os.path.basename(fastq) for fastq in fastqs]
    fastqs = [fastq.rsplit(".gz", 1)[0] for fastq in fastqs]
    fastqs = [fastq.rsplit(".fastq", 1)[0] for fastq in fastqs]

    final_output = expand(
        "results/qc/fastqc/{fastq_file}_fastqc.html",
        fastq_file=fastqs,
    )
    final_output.append("results/reports/multiqc/fastq.html")
    final_output.append(expand(
        "results/qc/samtools/idxstats/{sample_id}",
        sample_id=fastq_files["sample_id"].unique().tolist(),
    ))
    return final_output


def get_fastqc_input(wildcards):
    u = fastq_files.dropna()
    u["index"] = range(u.shape[0])
    u = pd.wide_to_long(u, ["fastq"], i="index", j="read")
    u["basename_prefix"] = [os.path.basename(x).rstrip("fastq.gz") for x in u.fastq]
    u = u.set_index(["basename_prefix"], drop=False)
    input_file = u.loc[wildcards.fastq_prefix, ["fastq"]].tolist()
    return input_file


def get_fastqc_reports():
    u = fastq_files.dropna()
    u["index"] = range(u.shape[0])
    u = pd.wide_to_long(u, ["fastq"], i="index", j="read")
    fastq_prefix = [os.path.basename(x).rstrip("fastq.gz") for x in u.fastq]
    fastq_reports = [
        f"results/qc/fastqc/{prefix}_fastqc.html" for prefix in fastq_prefix
    ]
    return fastq_reports


def get_hisat2_input(wildcards):
    u = fastq_files.dropna()
    fastq1 = fastq_files["fastq1"].groupby(
        fastq_files['sample_id']
    ).aggregate(
        lambda x: x.tolist()
    ).to_dict()
    fastq2 = fastq_files["fastq2"].groupby(
        fastq_files['sample_id']
    ).aggregate(
        lambda x: x.tolist()
    ).to_dict()
    return({
        "fastq1": fastq1[wildcards.sample_id],
        "fastq2": fastq2[wildcards.sample_id],
    })
