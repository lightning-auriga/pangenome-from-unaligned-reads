rule collect_insert_size_metrics:
    """
    Take initial alignment data for a library and compute
    approximate insert size distribution metrics.
    """
    input:
        bam="results/initial-alignments-sorted/{reference_genome}/all-pairs/{sampleid}.bam",
        bai="results/initial-alignments-sorted/{reference_genome}/all-pairs/{sampleid}.bai",
    output:
        histogram="results/insert-size-metrics/{reference_genome}/{sampleid}.insert_size_metrics.histogram.pdf",
        metrics="results/insert-size-metrics/{reference_genome}/{sampleid}.insert_size_metrics.tsv",
    params:
        java_args=config_resources["gatk-collectinsertsizemetrics"]["java-args"],
    benchmark:
        "results/performance_benchmarks/collect_insert_size_metrics/{reference_genome}/{sampleid}.tsv"
    conda:
        "../envs/gatk4.yaml"
    threads: config_resources["gatk-collectinsertsizemetrics"]["threads"]
    resources:
        mem_mb=config_resources["gatk-collectinsertsizemetrics"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["gatk-collectinsertsizemetrics"]["partition"],
            config_resources["partitions"],
        ),
        tmpdir=tempDir,
    shell:
        'gatk --java-options "{params.java_args}" CollectInsertSizeMetrics '
        "-INPUT {input.bam} -OUTPUT {output.metrics} -Histogram_FILE {output.histogram} "
        "-ASSUME_SORTED true -STOP_AFTER 2000000000"


rule rename_filter_fastq:
    """
    Emulate the same function from the masurca/superreads repo.
    Not sure why there are C++ and perl versions of this upstream.
    Kinda chaotic.

    Pulling in the existing samtools env just for bgzip.
    """
    input:
        R1="{prefix}_R1.fastq.gz",
        R2="{prefix}_R2.fastq.gz",
    output:
        fq="{prefix}.renamed.fastq.gz",
    benchmark:
        "results/performance_benchmarks/rename_filter_fastq/{prefix}.tsv"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    script:
        "../scripts/rename_filter_fastq.py"
