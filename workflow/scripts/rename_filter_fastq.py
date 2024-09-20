import gzip as gz
import re


def prepare_fastq_entry(f, library_tag, entry_count) -> str:
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
            val = "@{}{}\n".format(library_tag, entry_count)
        res += val[0:250]
    return res


def rename_filter_fastq(R1: str, R2: str, library_tag: str, combined: str) -> None:
    """
    Interleave paired reads from two fastqs into an output fastq.
    Not sure why masurca is doing it this way but maybe we'll find out..
    """
    r1_open = gz.open if R1.endswith(".gz") else open
    use_r2 = R2 is not None
    if use_r2:
        r2_open = gz.open if R2.endswith(".gz") else open
    out_open = gz.open if combined.endswith(".gz") else open
    with r1_open(R1, "rt") as f1:
        if use_r2:
            f2 = r2_open(R2, "rt")
        with out_open(combined, "wt") as out:
            entry_count = 1
            while 1:
                r1_output = prepare_fastq_entry(f1, library_tag, entry_count)
                if use_r2:
                    r2_output = prepare_fastq_entry(f2, library_tag, entry_count + 1)
                if use_r2 and (len(r1_output) > 0) != (len(r2_output) > 0):
                    raise ValueError("paired fastqs are of uneven lengths")
                if len(r1_output) == 0:
                    break
                out.write(r1_output)
                if use_r2:
                    out.write(r2_output)
                entry_count += 2 if use_r2 else 1
        if use_r2:
            f2.close()


rename_filter_fastq(
    snakemake.input["R1"],  # noqa: F821
    snakemake.input["R2"] if "R2" in snakemake.input.keys() else None,  # noqa: F821
    snakemake.params["library_tag"],  # noqa: F821
    snakemake.output["fq"],  # noqa: F821
)
