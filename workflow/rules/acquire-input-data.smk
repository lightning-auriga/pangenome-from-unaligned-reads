rule get_reference:
    """
    Get the configured reference genome in fasta format.
    """
    output:
        fn="results/references/{}/ref.fasta".format(config["reference-genome"]["name"]),
        tmp=temp("results/references/{}/ref.fasta.tmp".format(config["reference-genome"]["name"])),
    benchmark:
        "results/performance_benchmarks/get_reference/result.tsv"
    params:
        config["reference-genome"]["fasta"],
    conda:
        "../envs/awscli.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        'if [[ "{params}" == "s3://"* ]] ; then aws s3 cp {params} {output.tmp} ; '
        'elif [[ "{params}" == "http://"* ]] || [[ "{params}" == "https://"* ]] || [[ "{params}" == "ftp://"* ]] ; then wget -O {output.tmp} {params} ; '
        "else cp {params} {output.tmp} ; fi ; "
        'if [[ "{params}" = *".gz" ]] && [[ "{output.fn}" != *".gz" ]] ; then cat {output.tmp} | gunzip -c > {output.fn} ; '
        "else cp {output.tmp} {output.fn} ; fi"


rule samtools_index_reference:
    """
    Use samtools faidx to create fai index from reference fasta.
    """
    input:
        fasta="results/references/{reference_genome}/ref.fasta",
    output:
        fai="results/references/{reference_genome}/ref.fasta.fai",
    benchmark:
        "results/performance_benchmarks/samtools_index_reference/{reference_genome}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools faidx {input.fasta}"


rule bowtie2_index_reference:
    """
    Generate index for bowtie2 alignment from a reference fasta.
    """
    input:
        fasta="results/references/{reference_genome}/ref.fasta",
        fai="results/references/{reference_genome}/ref.fasta.fai",
    output:
        fwd_1="results/bowtie2-index/{reference_genome}/ref.fasta.1.bt2",
        fwd_2="results/bowtie2-index/{reference_genome}/ref.fasta.2.bt2",
        fwd_3="results/bowtie2-index/{reference_genome}/ref.fasta.3.bt2",
        fwd_4="results/bowtie2-index/{reference_genome}/ref.fasta.4.bt2",
        rev_1="results/bowtie2-index/{reference_genome}/ref.fasta.rev.1.bt2",
        rev_2="results/bowtie2-index/{reference_genome}/ref.fasta.rev.2.bt2",
    benchmark:
        "results/performance_benchmarks/bowtie2_index_reference/{reference_genome}.tsv"
    params:
        index_prefix="results/bowtie2-index/{reference_genome}/ref.fasta",
    conda:
        "../envs/bowtie2.yaml"
    threads: config_resources["bowtie2-build"]["threads"]
    resources:
        mem_mb=config_resources["bowtie2-build"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["bowtie2-build"]["partition"], config_resources["partitions"]
        ),
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.index_prefix}"


rule get_fastq:
    """
    Acquire a copy of a fastq specified in the manifest.
    """
    output:
        "results/input-fastqs/{sampleid}_{readgroup}.fastq.gz",
    benchmark:
        "results/performance_benchmarks/get_fastq/{sampleid}_{readgroup}.tsv"
    params:
        lambda wildcards: manifest.loc[wildcards.sampleid, wildcards.readgroup],
    conda:
        "../envs/awscli.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        'if [[ "{params}" == "s3://"* ]] ; then aws s3 cp {params} {output} ; '
        'elif [[ "{params}" == "http://"* ]] || [[ "{params}" == "https://"* ]] || [[ "{params}" == "ftp://"* ]] ; then wget -O {output} {params} ; '
        "else cp {params} {output} ; fi"
