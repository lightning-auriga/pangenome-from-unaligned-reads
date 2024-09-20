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


rule rename_filter_fastq_unmapped_reads:
    """
    Emulate the same function from the masurca/superreads repo.
    Not sure why there are C++ and perl versions of this upstream.
    Kinda chaotic.

    Pulling in the existing samtools env just for bgzip.
    """
    input:
        R1="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}_R1.fastq.gz",
        R2="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}_R2.fastq.gz",
    output:
        fq="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}.renamed.fastq.gz",
    params:
        library_tag="pe",
    benchmark:
        "results/performance_benchmarks/rename_filter_fastq_unmapped_reads/{reference_genome}/{sampleid}.tsv"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    script:
        "../scripts/rename_filter_fastq.py"


rule rename_filter_fastq_half_aligned_reads:
    """
    Emulate the same function from the masurca/superreads repo.
    Not sure why there are C++ and perl versions of this upstream.
    Kinda chaotic.

    Pulling in the existing samtools env just for bgzip.
    """
    input:
        R1="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}_unmapped.fastq.gz",
    output:
        fq="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}_unmapped.renamed.fastq.gz",
    params:
        library_tag="se",
    benchmark:
        "results/performance_benchmarks/rename_filter_fastq_half_aligned_reads/{reference_genome}/{sampleid}.tsv"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    script:
        "../scripts/rename_filter_fastq.py"


rule compute_read_metrics:
    """
    Take the upstream masurca logic for these calculations
    and spin it out into a python script for better handling
    of corner cases and unit testing.
    """
    input:
        "results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}.renamed.fastq.gz",
        "results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}_unmapped.renamed.fastq.gz",
    output:
        "results/read-metrics/{reference_genome}/{sampleid}.read_metrics.tsv",
    params:
        sampled_read_count=40000,
    benchmark:
        "results/performance_benchmarks/compute_read_metrics/{reference_genome}/{sampleid}.tsv"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    script:
        "../scripts/compute_read_metrics.py"
