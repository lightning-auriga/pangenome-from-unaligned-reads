import gzip as gz


def extract_read(f) -> str:
    """
    Extract a single fastq read (but devour four lines) from an input fastq file.
    """
    res = ""
    for i in range(4):
        val = f.readline()
        if val == "" and i != 0:
            raise ValueError("input fastq stream has improper line count")
        elif val == "":
            return ""
        if i == 1:
            res = val
    return res


def compute_average_read_length_from_file(fn: str, read_count: int) -> (int, int):
    """
    Compute numerator and denominator of read length average calculation from
    a single fastq file, filtering out reads of less than 31 bases.

    The assorted magic numbers in this script are set upstream
    without justification.
    """
    fq_open = gz.open if fn.endswith(".gz") else open
    with fq_open(fn, "rt") as f:
        num = 0
        denom = 0
        while 1:
            single_read = extract_read(f)
            if len(single_read) == 0:
                break
            if len(single_read) > 31:
                num += len(single_read) - 1
                denom += 1
            if denom >= read_count:
                break
        if denom < 1:
            raise ValueError(
                "insufficiently many valid reads available in fq {}: {}".format(fn, denom)
            )
        return (num, denom)


def compute_average_read_length_from_files(fns: list, read_count: int, out_fn: str) -> None:
    """
    Compute average read length from a set of fastq files by sampling the top
    N reads from each of the files and taking the average from those reads.

    This logic changes a few aspects of upstream to guarantee that a minimum
    number of reads is actually processed.

    Kinda wish I had boost accumulators right now.
    """
    total_num = 0
    total_denom = 0
    for fn in fns:
        num, denom = compute_average_read_length_from_file(fn, read_count)
        total_num += num
        total_denom += denom
    if total_denom < 1:
        raise ValueError("impossible denominator, cannot compute average read length")
    with open(out_fn, "w") as f:
        f.write("{}\t{}\t{}\n".format(total_num, total_denom, total_num / total_denom))


compute_average_read_length_from_files(
    snakemake.input,  # noqa: F821
    snakemake.params["sampled_read_count"],  # noqa: F821
    snakemake.output[0],  # noqa: F821
)
