from uuid import uuid4
from textwrap import dedent
from io import StringIO
import pandas as pd
import pytest
import sys

### VIASH START
meta = {
    "resources_dir": "./src/stats/generate_pool_statistics/",
    "executable": "target/executable/stats/generate_pool_statistics/generate_pool_statistics",
    "config": "src/stats/generate_pool_statistics/config.vsh.yaml"
}
### VIASH END

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper


@pytest.fixture
def random_tsv_path(random_path):
    def wrapper():
        return random_path(".tsv")
    return wrapper


@pytest.fixture
def simple_input_file_one(random_tsv_path, request):
    prefix = request.param
    mito_name = f"{prefix}M{'T' if not prefix else ''}"

    contents = dedent(
    f"""\
    WellBC	Chr	NumberOfReads	NumberOfGenes
    AGG	{prefix}1	2	1
    AGG	{prefix}2	3	2
    AGG	{prefix}3	4	2
    AGG	{mito_name}	4	2
    AGG	{prefix}X	2	3
    AGG	ERCC-1	1	1
    AGG	ERCC-2	1	1
    """)
    output_file = random_tsv_path()
    with output_file.open("w") as open_file:
        open_file.write(contents)
    return output_file


@pytest.fixture
def simple_input_file_two(random_tsv_path, request):
    prefix = request.param
    contents = dedent(
    f"""\
    WellBC	Chr	NumberOfReads	NumberOfGenes
    CCC	{prefix}2	2	1
    CCC	{prefix}3	3	2
    CCC	{prefix}5	4	2
    CCC	{prefix}1	4	2
    CCC	{prefix}Y	2	3
    CCC	{prefix}X	2	3
    CCC	ERCC-3	1	1
    CCC	ERCC-2	1	1
    """)
    output_file = random_tsv_path()
    with output_file.open("w") as open_file:
        open_file.write(contents)
    return output_file

@pytest.mark.parametrize("simple_input_file_one,simple_input_file_two,expected", [("chr", "chr", "chr"), ("", "", "")], 
                         indirect=["simple_input_file_one", "simple_input_file_two"])
def test_generate_pool_statistics_simple(run_component, simple_input_file_one,
                                         simple_input_file_two, random_tsv_path, expected):
    
    output_path = random_tsv_path()
    run_component([
        "--nrReadsNrGenesPerChrom", simple_input_file_one,
        "--nrReadsNrGenesPerChrom", simple_input_file_two,
        "--nrReadsNrGenesPerChromPool", output_path
    ])
    mito_name = f"{expected}M{'T' if not expected else ''}"
    expected_output = StringIO(dedent(
    f"""\
    WellBC	ERCC-1	ERCC-2	ERCC-3	{expected}1	{expected}2	{expected}3	{expected}5	{mito_name}	{expected}X	{expected}Y	SumReads	pctMT	pctERCC	pctChrom	NumberOfGenes
    AGG	1	1	0	2	3	4	0	4	2	0	17	23.53	11.76	52.94	12
    CCC	0	1	1	4	2	3	4	0	2	2	19	0.0	10.53	68.42	15
    """))
    assert output_path.is_file()
    contents = pd.read_csv(output_path, sep="\t")
    expected_frame = pd.read_csv(expected_output, sep="\t")
    pd.testing.assert_frame_equal(contents, expected_frame, check_like=True)


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))