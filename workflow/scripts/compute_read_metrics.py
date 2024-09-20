import gzip as gz

import numpy as np


def extract_read(f) -> (str, int):
    """
    Extract a single fastq read (but devour four lines) from an input fastq file.
    While thus computing, get the minimum qual score for the read, in an ad hoc
    effort to estimate which phred scale is in use.
    """
    res = ""
    min_qual = 128
    for i in range(4):
        val = f.readline()
        if val == "" and i != 0:
            raise ValueError("input fastq stream has improper line count")
        elif val == "":
            return ("", min_qual)
        if i == 1:
            res = val
        elif i == 3:
            min_qual = min([ord(x) for x in val])
    return (res, min_qual)


def compute_read_lengths_and_gc_from_file(fn: str, read_count: int) -> (int, int):
    """
    Compute various read metrics from a single fastq file.

    The assorted magic numbers in this script are set upstream
    without justification.
    """
    fq_open = gz.open if fn.endswith(".gz") else open
    read_lengths = []
    read_gc = []
    min_qual_ord = 128
    with fq_open(fn, "rt") as f:
        while 1:
            single_read, min_qual = extract_read(f)
            if len(single_read) == 0:
                break
            read_lengths.append(len(single_read) - 1)
            read_gc.append(len([1 for x in single_read if x in ["G", "C", "g", "c"]]))
            if min_qual < min_qual_ord:
                min_qual_ord = min_qual
            if len(read_lengths) >= read_count:
                break
    return (read_lengths, read_gc, min_qual_ord)


def compute_read_metrics_from_files(fns: list, read_count: int, out_fn: str) -> None:
    """
    Compute average read length, 25% read length, phred scale, GC content, and kmer length
    from a set of fastq files by sampling the top N reads from each of the files and
    estimating from those reads.

    This logic changes a few aspects of upstream to guarantee that a minimum
    number of reads is actually processed.
    """
    total_read_lengths = []
    total_read_gc = []
    min_qual_ord = 128
    for fn in fns:
        read_lengths, read_gc, min_qual = compute_read_lengths_and_gc_from_file(fn, read_count)
        total_read_lengths.extend(read_lengths)
        total_read_gc.extend(read_gc)
        if min_qual < min_qual_ord:
            min_qual_ord = min_qual
    frac_gc = sum(total_read_gc) / sum(total_read_lengths)
    av_len = np.mean([x for x in read_lengths if x > 31])
    min_len = np.quantile(read_lengths, 0.25)
    kmer = int(min_len * 0.66) if (frac_gc >= 0.35 and frac_gc <= 0.6) else int(min_len * 0.33)
    kmer += 1
    if kmer % 2 == 0:
        kmer += 1
    if kmer < 31:
        kmer = 31
    if kmer > 99:
        kmer = 99
    with open(out_fn, "w") as f:
        f.write(
            "{}\t{}\t{}\t{}\t{}\n".format(
                av_len, min_len, frac_gc, kmer, 33 if min_qual_ord < 64 else 64
            )
        )


compute_read_metrics_from_files(
    snakemake.input,  # noqa: F821
    snakemake.params["sampled_read_count"],  # noqa: F821
    snakemake.output[0],  # noqa: F821
)
