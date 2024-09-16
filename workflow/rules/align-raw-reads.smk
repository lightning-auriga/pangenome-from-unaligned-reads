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
        bam="results/initial-alignments/{reference_genome}/all-pairs/{sampleid}.bam",
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


rule extract_unaligned_pairs:
    """
    Take initial basic bowtie2 alignments and select read pairs that had neither read map
    to the corresponding reference genome.
    """
    input:
        bam="results/initial-alignments/{reference_genome}/all-pairs/{sampleid}.bam",
    output:
        bam="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}.bam",
    benchmark:
        "results/performance_benchmarks/extract_unaligned_pairs/{reference_genome}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools view --no-PG -e 'flag.unmap && flag.munmap' -o {output.bam} {input.bam}"


rule extract_half_aligned_pairs:
    """
    Take initial basic bowtie2 alignments and select read pairs that had exactly one read
    of the pair map to the corresponding reference genome.
    """
    input:
        bam="results/initial-alignments/{reference_genome}/all-pairs/{sampleid}.bam",
    output:
        bam="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}.bam",
    benchmark:
        "results/performance_benchmarks/extract_half_aligned_pairs/{reference_genome}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools view --no-PG -e '(flag.unmap && !flag.munmap) || (!flag.unmap && flag.munmap)' -o {output.bam} {input.bam}"


rule unaligned_pairs_bam_to_fastq:
    input:
        bam="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}.bam",
    output:
        fastq_r1="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}_R1.fastq.gz",
        fastq_r2="results/initial-alignments/{reference_genome}/unaligned-pairs/{sampleid}_R2.fastq.gz",
    benchmark:
        "results/performance_benchmarks/unaligned_pairs_bam_to_fastq/{reference_genome}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools fastq -1 {output.fastq_r1} -2 {output.fastq_r2} -0 /dev/null -s /dev/null {input.bam}"


rule half_aligned_pairs_bam_to_unmapped_read_fastq:
    input:
        bam="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}.bam",
    output:
        fastq="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}_unmapped.fastq.gz",
    benchmark:
        "results/performance_benchmarks/half_aligned_pairs_bam_to_unmapped_read_fastq/{reference_genome}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools view -e 'flag.unmap' -b {input.bam} | "
        "samtools fastq - | "
        "bgzip -c > {output.fastq}"


rule half_aligned_pairs_bam_to_mapped_read_fastq:
    input:
        bam="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}.bam",
    output:
        fastq="results/initial-alignments/{reference_genome}/half-aligned-pairs/{sampleid}_mapped.fastq.gz",
    benchmark:
        "results/performance_benchmarks/half_aligned_pairs_bam_to_mapped_read_fastq/{reference_genome}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools view -e '!flag.unmap' -b {input.bam} | "
        "samtools fastq - | "
        "bgzip -c > {output.fastq}"


rule sort_aligned_raw_reads:
    """
    As part of testing out interactions with the primary alignments, this rule sorts the alignments
    to allow indexing and random region access.
    """
    input:
        bam="results/initial-alignments/{reference_genome}/{pairset}/sampleid}.bam",
    output:
        bam="results/initial-alignments-sorted/{reference_genome}/{pairset}/{sampleid}.bam",
    benchmark:
        "results/performance_benchmarks/sort_aligned_raw_reads/{reference_genome}/{pairset}/{sampleid}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["samtools-sort"]["threads"]
    resources:
        mem_mb=config_resources["samtools-sort"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["samtools-sort"]["partition"], config_resources["partitions"]
        ),
        tmpdir=tempDir,
    shell:
        "samtools sort -@ {threads} -T {resources.tmpdir} -o {output.bam} {input.bam}"


rule samtools_index_bam:
    """
    Take an arbitrary bam file, which implicitly must be sorted, and index it.
    """
    input:
        bam="{prefix}.bam",
    output:
        bai="{prefix}.bai",
    benchmark:
        "results/performance_benchmarks/samtools_index_bam/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools index -@ {threads} -o {output.bai} {input.bam}"
