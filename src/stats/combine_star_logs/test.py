import pytest
import sys
import re
import pandas as pd
from pathlib import Path
from uuid import uuid4
from subprocess import CalledProcessError

### VIASH START
meta = {
    "resources_dir": "./src/stats/combine_star_logs/",
    "executable": "target/executable/stats/combine_star_logs/combine_star_logs",
    "config": "src/stats/combine_star_logs/config.vsh.yaml"
}
### VIASH END

@pytest.fixture
def test_resources_path():
    return Path(meta["resources_dir"]) / "test_data"

@pytest.fixture
def barcode_1_star_log(test_resources_path):
    return test_resources_path / "barcode_1" / "Log.final.out"

@pytest.fixture
def barcode_1_reads_per_gene_file(test_resources_path):
    return test_resources_path / "barcode_1" / "ReadsPerGene.out.tab"

@pytest.fixture
def barcode_1_summary(test_resources_path):
    return test_resources_path / "barcode_1" / "summary.csv"

@pytest.fixture
def barcode_2_star_log(test_resources_path):
    return test_resources_path / "barcode_2" / "Log.final.out"

@pytest.fixture
def barcode_2_reads_per_gene_file(test_resources_path):
    return test_resources_path / "barcode_2" / "ReadsPerGene.out.tab"

@pytest.fixture
def barcode_2_summary(test_resources_path):
    return test_resources_path / "barcode_2" / "summary.csv"

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper 

def test_incorrect_number_of_inputs_raises(run_component,
                                           barcode_1_star_log, barcode_2_star_log,
                                           barcode_1_reads_per_gene_file, barcode_2_reads_per_gene_file,
                                           barcode_1_summary, barcode_2_summary,
                                           random_path):
    output_path = random_path("txt")
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--barcodes", "foo;bar",
            "--star_logs", f"{barcode_1_star_log}", 
            "--reads_per_gene_logs", f"{barcode_1_reads_per_gene_file};{barcode_2_reads_per_gene_file}",
            "--gene_summary_logs", f"{barcode_1_summary};{barcode_2_summary}",
            "--output", output_path,
        ])
    assert re.search(r"ValueError: Expected the same number of inputs for 'star_logs' \(1\), "
                     r"'gene_summary_logs' \(2\), 'reads_per_gene_logs' \(2\) and 'barcodes' \(2\)\.",
            err.value.stdout.decode('utf-8'))



def test_equal_number_of_argument(run_component,
                                  barcode_1_star_log, barcode_2_star_log,
                                  barcode_1_reads_per_gene_file, barcode_2_reads_per_gene_file,
                                  barcode_1_summary, barcode_2_summary,
                                  random_path):
    output_path = random_path("txt")
    run_component([
        "--barcodes", "foo;bar",
        "--star_logs", f"{barcode_1_star_log};{barcode_2_star_log}", 
        "--reads_per_gene_logs", f"{barcode_1_reads_per_gene_file};{barcode_2_reads_per_gene_file}",
        "--gene_summary_logs", f"{barcode_1_summary};{barcode_2_summary}",
        "--output", output_path,
    ])
    expected_dict = {
        'NumberOfInputReads': [96398.0, 10155.0], 
        'NumberOfMappedReads': [70824.0, 7179.0], 
        'PctMappedReads': [73.47, 70.69], 
        'NumberOfReadsMappedToMultipleLoci': [0.0, 0.0], 
        'PectOfReadsMappedToMultipleLoci': [0.0, 0.0], 
        'NumberOfReadsMappedToTooManyLoci': [22281.0, 2248.0],
        'PectOfReadsMappedToTooManyLoci': [23.11, 22.14],
        'NumberOfReadsUnmappedTooManyMismatches': [0.0, 0.0], 
        'PectOfReadsUnmappedTooManyMismatches': [0.0, 0.0], 
        'NumberOfReadsUnmappedTooShort': [2697.0, 553.0], 
        'PectOfReadsUnmappedTooShort': [2.8, 5.45], 
        'NumberOfReadsUnmappedOther': [596.0, 175.0], 
        'PectOfReadsUnmappedOther': [0.62, 1.72], 
        'ReadsWithValidBarcodes': [0.999782, 0.999803],
        'SequencingSaturation': [0.0602963, 0.0539344], 
        'EstimatedNumberOfCells': [1.0, 1.0], 
        'NumberOfUMIs': [50370.0, 4701.0], 
        'NumberOfGenes': [8767.0, 2397.0],
        'NumberOfCountedReads': [17, 15],
    }
    expected = pd.DataFrame.from_dict(expected_dict)
    expected.index = pd.Index(["foo", "bar"], name="WellBC")
    assert output_path.is_file()

    contents = pd.read_csv(output_path, sep="\t", index_col=0)
    assert set(("NumberOfInputReads", "SequencingSaturation",
                "NumberOfGenes", "NumberOfUMIs", "NumberOfCountedReads",
                "PctMappedReads")).issubset(set(contents.columns))
    pd.testing.assert_frame_equal(contents, expected)



if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))