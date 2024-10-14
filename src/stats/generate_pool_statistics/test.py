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
    expected_dict = {
        "WellBC": ["AGG", "CCC"],
        "ERCC-1": ["1", "0"],
        "ERCC-2": ["1", "1"],
        "ERCC-3": ["0", "1"],
        f"{expected}1": ["2", "4"],
        f"{expected}2": ["3", "2"],
        f"{expected}3": ["4", "3"],
        f"{expected}5": ["0", "4"],
        f"{mito_name}": ["4", "0"],
        f"{expected}X": ["2", "2"],
        f"{expected}Y": ["0", "2"],
        "SumReads": ["17", "19"],
        "pctMT": ["23.53", "0"],
        "pctERCC": ["11.76", "10.53"],
        "pctChrom": ["52.94", "68.42"],
        "NumberOfGenes": ["12", "15"],
        "NumberOfMTReads": ["4", "0"],
        "NumberOfChromReads": ["9", "13"],
        "NumberOfERCCReads": ["2", "2"],
    }
    expected_frame = pd.DataFrame.from_dict(expected_dict, dtype=pd.StringDtype())
    assert output_path.is_file()
    contents = pd.read_csv(output_path, sep="\t", dtype=pd.StringDtype())
    pd.testing.assert_frame_equal(contents, expected_frame, check_like=True)


def test_only_numerical_chromosomes(run_component, random_tsv_path):
    """
    The chromosome column might be read as an integer instead of a string,
    make sure that a numerical column only works.
    """
    output_path = random_tsv_path()
    contents1 = dedent(
    f"""\
    WellBC	Chr	NumberOfReads	NumberOfGenes
    CCC	2	2	1
    CCC	3	3	2
    CCC	5	4	2
    CCC	1	4	2
    """)
    input_file_1 = random_tsv_path()
    with input_file_1.open("w") as open_file:
        open_file.write(contents1)

    contents2 = dedent(
    f"""\
    WellBC	Chr	NumberOfReads	NumberOfGenes
    AGG	2	2	1
    AGG	3	3	2
    AGG	5	4	2
    AGG	1	4	2
    """)
    input_file_2 = random_tsv_path()
    with input_file_2.open("w") as open_file:
        open_file.write(contents2)
        output_path = random_tsv_path()
    run_component([
        "--nrReadsNrGenesPerChrom", input_file_1,
        "--nrReadsNrGenesPerChrom", input_file_2,
        "--nrReadsNrGenesPerChromPool", output_path
    ])

    expected_dict = {
        "WellBC": ["AGG", "CCC"],
        "1": ["4", "4"],
        "2": ["2", "2"],
        "3": ["3", "3"],
        "5": ["4", "4"],
        "pctChrom": ["100", "100"],
        "pctMT": ["0", "0"],
        "pctERCC": ["0", "0"],
        "SumReads": ["13", "13"],
        "NumberOfGenes": ["7", "7"],
        "NumberOfERCCReads": ["0", "0"],
        "NumberOfChromReads": ["13", "13"],
        "NumberOfMTReads": ["0", "0"],
    }
    expected_frame = pd.DataFrame.from_dict(expected_dict,
                                            dtype=pd.StringDtype())

    assert output_path.is_file()
    contents = pd.read_csv(output_path, sep="\t", dtype=pd.StringDtype())
    pd.testing.assert_frame_equal(contents, expected_frame, check_like=True)


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))