import gzip as gz
import re


def prepare_fastq_entry(f, entry_count) -> str:
    """
    Extract a single fastq set (four lines) from an input fastq file.
    """
    res = ""
    for i in range(4):
        val = f.readline()
        if val == "" and i != 0:
            raise ValueError("input fastq stream has improper line count")
        elif val == "":
            return ""
        if i == 0:
            val = "@pe{}\n".format(entry_count)
        res += val
    return res


def rename_filter_fastq(R1: str, R2: str, combined: str) -> None:
    """
    Interleave paired reads from two fastqs into an output fastq.
    Not sure why masurca is doing it this way but maybe we'll find out..
    """
    r1_open = gz.open if R1.endswith(".gz") else open
    r2_open = gz.open if R2.endswith(".gz") else open
    out_open = gz.open if combined.endswith(".gz") else open
    with r1_open(R1, "rt") as f1:
        with r2_open(R2, "rt") as f2:
            with out_open(combined, "wt") as out:
                entry_count = 1
                while 1:
                    r1_output = prepare_fastq_entry(f1, entry_count)
                    r2_output = prepare_fastq_entry(f2, entry_count + 1)
                    if (len(r1_output) > 0) != (len(r2_output) > 0):
                        raise ValueError("paired fastqs are of uneven lengths")
                    if len(r1_output) == 0:
                        break
                    out.write(r1_output)
                    out.write(r2_output)
                    entry_count += 2


rename_filter_fastq(
    snakemake.input["R1"],  # noqa: F821
    snakemake.input["R2"],  # noqa: F821
    snakemake.output["fq"],  # noqa: F821
)
