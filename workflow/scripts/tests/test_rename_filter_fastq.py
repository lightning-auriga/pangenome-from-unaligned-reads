import gzip as gz
import os
import pathlib
import re
import runpy

import pandas as pd
import pytest
import snakemake.script as sms
import yaml
from snakemake.io import Namedlist


@pytest.fixture
def r1_fastq_filename(tmp_path):
    return tmp_path / "rename_filter_fastq" / "paired-end-R1.fastq.gz"


@pytest.fixture
def r2_fastq_filename(tmp_path):
    return tmp_path / "rename_filter_fastq" / "paired-end-R2.fastq.gz"


@pytest.fixture
def single_fastq_filename(tmp_path):
    return tmp_path / "rename_filter_fastq" / "single-end.fastq.gz"


@pytest.fixture
def output_paired_fastq_filename(tmp_path):
    return tmp_path / "rename_filter_fastq" / "output" / "paired-end.fastq.gz"


@pytest.fixture
def output_single_fastq_filename(tmp_path):
    return tmp_path / "rename_filter_fastq" / "output" / "single-end.fastq.gz"


@pytest.fixture
def paired_end_r1_content():
    res = ["@R01", "ACTT", "+", "FFFF", "@R02", "TTAA", "+", "FEAC", "@R03", "CCCC", "+", "####"]
    return "\n".join(res) + "\n"


@pytest.fixture
def paired_end_r2_content():
    res = ["@R04", "CAAC", "+", "DDDD", "@R05", "TAAA", "+", "CCCC", "@R06", "AAAA", "+", "BBBC"]
    return "\n".join(res) + "\n"


@pytest.fixture
def single_end_content():
    res = ["@R07", "AACC", "+", "AACC", "@R08", "CACA", "+", "BAFE"]
    return "\n".join(res) + "\n"


@pytest.fixture
def paired_library_tag():
    return "pe"


@pytest.fixture
def single_library_tag():
    return "se"


@pytest.fixture
def expected_single_output(single_end_content, single_library_tag):
    return "@se1\nAACC\n+\nAACC\n@se2\nCACA\n+\nBAFE\n"


@pytest.fixture
def expected_paired_output(paired_end_r1_content, paired_end_r2_content, paired_library_tag):
    return (
        "@pe1\nACTT\n+\nFFFF\n@pe2\nCAAC\n+\nDDDD\n"
        + "@pe3\nTTAA\n+\nFEAC\n@pe4\nTAAA\n+\nCCCC\n"
        + "@pe5\nCCCC\n+\n####\n@pe6\nAAAA\n+\nBBBC\n"
    )


def write_fastq_file(fn, content):
    pathlib.Path(os.path.dirname(fn)).mkdir(parents=True, exist_ok=True)
    with gz.open(fn, "wt") as f:
        f.write(content)


@pytest.fixture
def paired_end_files(
    r1_fastq_filename, paired_end_r1_content, r2_fastq_filename, paired_end_r2_content
):
    write_fastq_file(r1_fastq_filename, paired_end_r1_content)
    write_fastq_file(r2_fastq_filename, paired_end_r2_content)
    return (r1_fastq_filename, r2_fastq_filename)


@pytest.fixture
def single_end_file(single_fastq_filename, single_end_content):
    write_fastq_file(single_fastq_filename, single_end_content)
    return single_fastq_filename


@pytest.fixture
def snakemake_input_paired_end(paired_end_files):
    res = Namedlist(fromdict={"R1": str(paired_end_files[0]), "R2": str(paired_end_files[1])})
    return res


@pytest.fixture
def snakemake_input_single_end(single_end_file):
    res = Namedlist(fromdict={"R1": str(single_end_file)})
    return res


@pytest.fixture
def snakemake_output_paired_end(output_paired_fastq_filename):
    res = Namedlist(fromdict={"fq": str(output_paired_fastq_filename)})
    return res


@pytest.fixture
def snakemake_output_single_end(output_single_fastq_filename):
    res = Namedlist(fromdict={"fq": str(output_single_fastq_filename)})
    return res


@pytest.fixture
def snakemake_params_paired_end(paired_library_tag):
    res = Namedlist(fromdict={"library_tag": paired_library_tag})
    return res


@pytest.fixture
def snakemake_params_single_end(single_library_tag):
    res = Namedlist(fromdict={"library_tag": single_library_tag})
    return res


@pytest.fixture
def snakemake_object_paired(
    snakemake_input_paired_end, snakemake_output_paired_end, snakemake_params_paired_end
):
    res = sms.Snakemake(
        snakemake_input_paired_end,
        snakemake_output_paired_end,
        snakemake_params_paired_end,
        Namedlist(),
        1,
        Namedlist(),
        Namedlist(),
        {},
        "",
        [],
    )
    res.input = snakemake_input_paired_end
    res.output = snakemake_output_paired_end
    res.params = snakemake_params_paired_end
    return res


@pytest.fixture
def snakemake_object_single(
    snakemake_input_single_end, snakemake_output_single_end, snakemake_params_single_end
):
    res = sms.Snakemake(
        snakemake_input_single_end,
        snakemake_output_single_end,
        snakemake_params_single_end,
        Namedlist(),
        1,
        Namedlist(),
        Namedlist(),
        {},
        "",
        [],
    )
    res.input = snakemake_input_single_end
    res.output = snakemake_output_single_end
    res.params = snakemake_params_single_end
    return res


def test_rename_filter_fastq_paired_end(snakemake_object_paired, expected_paired_output):
    pathlib.Path(os.path.dirname(snakemake_object_paired.output[0])).mkdir(
        parents=True, exist_ok=True
    )
    runpy.run_path(
        "workflow/scripts/rename_filter_fastq.py",
        init_globals={"snakemake": snakemake_object_paired},
    )
    assert pathlib.Path(snakemake_object_paired.output[0]).is_file()
    with gz.open(snakemake_object_paired.output[0], "rt") as f:
        observed = "".join(f.readlines())
    expected = expected_paired_output
    assert observed == expected


def test_rename_filter_fastq_single_end(snakemake_object_single, expected_single_output):
    pathlib.Path(os.path.dirname(snakemake_object_single.output[0])).mkdir(
        parents=True, exist_ok=True
    )
    runpy.run_path(
        "workflow/scripts/rename_filter_fastq.py",
        init_globals={"snakemake": snakemake_object_single},
    )
    assert pathlib.Path(snakemake_object_single.output[0]).is_file()
    with gz.open(snakemake_object_single.output[0], "rt") as f:
        observed = "".join(f.readlines())
    expected = expected_single_output
    assert observed == expected
