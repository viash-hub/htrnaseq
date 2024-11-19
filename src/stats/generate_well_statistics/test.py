import sys
import pytest
import pysam
from uuid import uuid4
from pathlib import Path
from textwrap import dedent

### VIASH START
meta = {
    "resources_dir": "./src/stats/generate_well_statistics/",
    "executable": "target/executable/stats/generate_well_statistics/generate_well_statistics",
    "config": "src/stats/generate_well_statistics/config.vsh.yaml"
}
### VIASH END

def assert_file_content_equals(file_to_check, expected):
    with file_to_check.open('r') as open_file:
        contents = open_file.read()
        assert contents == expected


@pytest.fixture
def input_sam_path():
    return Path(meta["resources_dir"]) / "test.sam"


@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper 

@pytest.fixture
def random_bam_path(random_path):
    def wrapper():
        return random_path(".bam")
    return wrapper


@pytest.fixture
def sam_to_bam(random_bam_path):
    def wrapper(sam_file):
        out_path = random_bam_path()
        with pysam.AlignmentFile(sam_file, "r") as infile, \
            pysam.AlignmentFile(out_path, "wb", template=infile) as outfile:
            for s in infile:
                outfile.write(s)
        infile.close()
        return out_path
    return wrapper


def test_generate_well_statistics_simple_bam(run_component, input_sam_path, sam_to_bam, random_path):
    bam_file = sam_to_bam(input_sam_path)
    processed_bam = random_path("tsv")
    reads_per_chromosome = random_path("tsv")
    nr_reads_nr_umis_per_cb = random_path("tsv")
    top_onehundred_umis = random_path("tsv")
    run_component([
        "--input", bam_file,
        "--processedBAMFile", processed_bam,
        "--nrReadsNrGenesPerChrom", reads_per_chromosome,
        "--nrReadsNrUMIsPerCB", nr_reads_nr_umis_per_cb,
        "--umiFreqTop", top_onehundred_umis,
        "--barcode", "ACGT",
        "--well_id", "A1",
    ])
    for file_path in (processed_bam, reads_per_chromosome,
                      nr_reads_nr_umis_per_cb, top_onehundred_umis):
        assert file_path.is_file()

    expected_processed_bam = \
    dedent("""\
    WellBC	WellID	Chr	CB	UB	GX	GN
    ACGT	A1	1	ACA	CGG	gene1	gene1
    ACGT	A1	1	ACA	CGG	gene1	gene1
    ACGT	A1	2	GGG	GTT	gene2	gene2
    ACGT	A1	2	GGG	GTC	gene3	gene3
    """)

    expected_reads_per_chromosome = \
    dedent("""\
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    ACGT	A1	1	2	1
    ACGT	A1	2	2	2
    """)

    expected_nr_reads_nr_umis_per_cb = \
    dedent("""\
    WellBC	WellID	CB	NumberOfReads	nrUMIs
    ACGT	A1	ACA	2	1
    ACGT	A1	GGG	2	2
    """)

    expected_top_onehundred_umis = \
    dedent("""\
    WellBC	WellID	UB	N
    ACGT	A1	CGG	2
    ACGT	A1	GTC	1
    ACGT	A1	GTT	1
    """)

    assert_file_content_equals(processed_bam, expected_processed_bam)
    assert_file_content_equals(reads_per_chromosome, expected_reads_per_chromosome)
    assert_file_content_equals(nr_reads_nr_umis_per_cb, expected_nr_reads_nr_umis_per_cb)
    assert_file_content_equals(top_onehundred_umis, expected_top_onehundred_umis)


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))