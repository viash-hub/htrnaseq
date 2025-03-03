from uuid import uuid4
from textwrap import dedent
from subprocess import CalledProcessError
import pandas as pd
import re
import pytest
import sys
from pathlib import Path

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
        return Path(tmp_path / f"{uuid4()}{extension}")
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
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    AGG	A1	{prefix}1	2	1
    AGG	A1	{prefix}2	3	2
    AGG	A1	{prefix}3	4	2
    AGG	A1	{mito_name}	4	2
    AGG	A1	{prefix}X	2	3
    AGG	A1	ERCC-1	1	1
    AGG	A1	ERCC-2	1	1
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
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    CCC	B2	{prefix}2	2	1
    CCC	B2	{prefix}3	3	2
    CCC	B2	{prefix}5	4	2
    CCC	B2	{prefix}1	4	2
    CCC	B2	{prefix}Y	2	3
    CCC	B2	{prefix}X	2	3
    CCC	B2	ERCC-3	1	1
    CCC	B2	ERCC-2	1	1
    """)
    output_file = random_tsv_path()
    with output_file.open("w") as open_file:
        open_file.write(contents)
    return output_file

@pytest.fixture
def empty_input_file(random_tsv_path):
    contents = dedent(
    f"""\
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
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
        "WellID": ["A1", "B2"],
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
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    CCC	B2	2	2	1
    CCC	B2	3	3	2
    CCC	B2	5	4	2
    CCC	B2	1	4	2
    """)
    input_file_1 = random_tsv_path()
    with input_file_1.open("w") as open_file:
        open_file.write(contents1)

    contents2 = dedent(
    f"""\
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    AGG	A1	2	2	1
    AGG	A1	3	3	2
    AGG	A1	5	4	2
    AGG	A1	1	4	2
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
        "WellID": ["A1", "B2"],
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


@pytest.mark.parametrize("simple_input_file_one", [("")],
                         indirect=["simple_input_file_one"])
def test_empty_input_raises(run_component, simple_input_file_one, empty_input_file, random_tsv_path):
    """
    When an input file contains no data, raise an error.
    """
    output_path = random_tsv_path()
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--nrReadsNrGenesPerChrom", simple_input_file_one,
            "--nrReadsNrGenesPerChrom", empty_input_file,
            "--nrReadsNrGenesPerChromPool", output_path
        ])
    assert re.search(
        rf"{empty_input_file.name} does not seem to contain any information",
        err.value.stdout.decode("utf-8"),
    )

def test_remove_chromosomes_with_no_counts(run_component, random_tsv_path):
    """
    If a chromosome has no counts across all of the wells, it should
    not be included in the output
    """
    output_path = random_tsv_path()
    contents1 = dedent(
    f"""\
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    CCC	B2	2	2	1
    CCC	B2	3	3	2
    CCC	B2	5	4	2
    CCC	B2	1	4	2
    CCC	B2	empty	0	0
    """)
    input_file_1 = random_tsv_path()
    with input_file_1.open("w") as open_file:
        open_file.write(contents1)

    contents2 = dedent(
    f"""\
    WellBC	WellID	Chr	NumberOfReads	NumberOfGenes
    AGG	A1	2	2	1
    AGG	A1	3	3	2
    AGG	A1	5	4	2
    AGG	A1	1	4	2
    AGG	A1	empty	0	0
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
    # Here, the chromosome called "empty" should not be included
    expected_dict = {
        "WellBC": ["AGG", "CCC"],
        "WellID": ["A1", "B2"],
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