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
