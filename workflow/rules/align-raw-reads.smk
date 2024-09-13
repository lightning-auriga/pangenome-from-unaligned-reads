rule bowtie2_align_raw_reads:
    """
    Take input fastq data and align to the standard reference
    genome. Based on the upstream reference, this alignment
    should be conducted with bowtie2.
    """
    input:
        R1="results/input-fastqs/{sampleid}_R1.fastq.gz",
        R2="results/input-fastqs/{sampleid}_R2.fastq.gz",
        fai="results/references/{reference_genome}/ref.fasta.fai",
        index=expand(
            "results/bowtie2-index/{{reference_genome}}/ref.fasta.{suffix}",
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"],
        ),
    output:
        bam="results/initial-alignments/{reference_genome}/{sampleid}.bam",
    benchmark:
        "results/performance_benchmarks/bowtie2_align_raw_reads/{reference_genome}/{sampleid}.tsv"
    params:
        index_prefix="results/bowtie2-index/{reference_genome}/ref.fasta",
    conda:
        "../envs/bowtie2.yaml"
    threads: config_resources["bowtie2-align"]["threads"]
    resources:
        mem_mb=config_resources["bowtie2-align"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["bowtie2-align"]["partition"], config_resources["partitions"]
        ),
    shell:
        "bowtie2 -x {params.index_prefix} -1 {input.R1} -2 {input.R2} --threads {threads} | "
        "samtools view -@1 -bt {input.fai} -o {output.bam} - "
