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
def barcode_1_features_stats(test_resources_path):
    return test_resources_path / "barcode_1" / "Features.stats"

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
def barcode_2_features_stats(test_resources_path):
    return test_resources_path / "barcode_2" / "Features.stats"

@pytest.fixture
def no_reads_mapped_star_log(test_resources_path):
    return test_resources_path / "empty" / "Log.final.out"

@pytest.fixture
def no_reads_mapped_reads_per_gene_file(test_resources_path):
    return test_resources_path / "empty" / "ReadsPerGene.out.tab"

@pytest.fixture
def no_reads_mapped_summary(test_resources_path):
    return test_resources_path / "empty" / "summary.csv"

@pytest.fixture
def no_reads_mapped_features_stats(test_resources_path):
    return test_resources_path / "empty" / "Features.stats"

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper 

def test_incorrect_number_of_inputs_raises(run_component,
                                           barcode_1_star_log,
                                           barcode_1_reads_per_gene_file, barcode_2_reads_per_gene_file,
                                           barcode_1_summary, barcode_2_summary,
                                           barcode_1_features_stats, barcode_2_features_stats,
                                           random_path):
    output_path = random_path("txt")
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--barcodes", "foo;bar",
            "--star_logs", f"{barcode_1_star_log}", 
            "--reads_per_gene_logs", f"{barcode_1_reads_per_gene_file};{barcode_2_reads_per_gene_file}",
            "--gene_summary_logs", f"{barcode_1_summary};{barcode_2_summary}",
            "--features_stats", f"{barcode_1_features_stats};{barcode_2_features_stats}",
            "--output", output_path,
        ])
    assert re.search(r"ValueError: Expected the same number of inputs for 'star_logs' \(1\), "
                     r"'gene_summary_logs' \(2\), 'reads_per_gene_logs' \(2\), 'features_stats' "
                     r"\(2\) and 'barcodes' \(2\)\.",
            err.value.stdout.decode('utf-8'))



def test_equal_number_of_argument(run_component,
                                  barcode_1_star_log, barcode_2_star_log,
                                  barcode_1_reads_per_gene_file, barcode_2_reads_per_gene_file,
                                  barcode_1_features_stats, barcode_2_features_stats,
                                  barcode_1_summary, barcode_2_summary,
                                  random_path):
    output_path = random_path("txt")
    run_component([
        "--barcodes", "foo;bar",
        "--star_logs", f"{barcode_1_star_log};{barcode_2_star_log}", 
        "--reads_per_gene_logs", f"{barcode_1_reads_per_gene_file};{barcode_2_reads_per_gene_file}",
        "--gene_summary_logs", f"{barcode_1_summary};{barcode_2_summary}",
        "--features_stats", f"{barcode_1_features_stats};{barcode_2_features_stats}",
        "--output", output_path,
    ])
    # We use strings here to make a comparison of the file contents without
    # doing any inferences of the numerical data type (i.e. exact file contents).
    expected_dict = {
        'NumberOfInputReads': ["96398", "10155"], 
        'NumberOfMappedReads': ["70824", "7179"], 
        'PctMappedReads': ["73.47", "70.69"], 
        'NumberOfReadsMappedToMultipleLoci': ["0", "0"], 
        'PectOfReadsMappedToMultipleLoci': ["0", "0"], 
        'NumberOfReadsMappedToTooManyLoci': ["22281", "2248"],
        'PectOfReadsMappedToTooManyLoci': ["23.11", "22.14"],
        'NumberOfReadsUnmappedTooManyMismatches': ["0", "0"], 
        'PectOfReadsUnmappedTooManyMismatches': ["0", "0"], 
        'NumberOfReadsUnmappedTooShort': ["2697", "553"], 
        'PectOfReadsUnmappedTooShort': ["2.8", "5.45"], 
        'NumberOfReadsUnmappedOther': ["596", "175"], 
        'PectOfReadsUnmappedOther': ["0.62", "1.72"], 
        'ReadsWithValidBarcodes': ["0.999782", "0.999803"],
        'SequencingSaturation': ["0.0602963", "0.0539344"], 
        'Q30BasesInCB+UMI': ["0.980096", "0.984461"],
        'ReadsMappedToTranscriptome:Unique+MultipeGenes': ["0.60411", "0.530871"],
        'EstimatedNumberOfCells': ["1", "1"],
        'FractionOfReadsInCells': ["1", "1"],
        'MeanReadsPerCell': ["53602", "4969"],
        'NumberOfUMIs': ["50370", "4701"], 
        'NumberOfGenes': ["8767", "2397"],
        'NumberOfCountedReads': ["17", "15"],
    }
    expected = pd.DataFrame.from_dict(expected_dict, dtype=pd.StringDtype())
    expected.index = pd.Index(["foo", "bar"], name="WellBC", dtype=pd.StringDtype())
    assert output_path.is_file()

    contents = pd.read_csv(output_path, sep="\t", index_col=0, dtype=pd.StringDtype())
    assert set(("NumberOfInputReads", "SequencingSaturation",
                "NumberOfGenes", "NumberOfUMIs", "NumberOfCountedReads",
                "PctMappedReads")).issubset(set(contents.columns))
    pd.testing.assert_frame_equal(contents, expected)

def test_empty(run_component, no_reads_mapped_star_log,
               no_reads_mapped_reads_per_gene_file, no_reads_mapped_summary,
               no_reads_mapped_features_stats,
               random_path):
    """
    Sometimes the summary.csv contains '-nan' values, make sure they
    are properly handled.
    """
    output_path = random_path("txt")
    run_component([
        "--barcodes", "foo",
        "--star_logs", no_reads_mapped_star_log,
        "--reads_per_gene_logs", no_reads_mapped_reads_per_gene_file,
        "--gene_summary_logs", no_reads_mapped_summary,
        "--features_stats", no_reads_mapped_features_stats,
        "--output", output_path,
    ])
    expected_dict = {
        'NumberOfInputReads': ["1327"],
        'NumberOfMappedReads': ["116"],
        'PctMappedReads': ["8.74"],
        'NumberOfReadsMappedToMultipleLoci': ["0"],
        'PectOfReadsMappedToMultipleLoci': ["0"],
        'NumberOfReadsMappedToTooManyLoci': ["43"],
        'PectOfReadsMappedToTooManyLoci': ["3.24"],
        'NumberOfReadsUnmappedTooManyMismatches': ["0"],
        'PectOfReadsUnmappedTooManyMismatches': ["0"],
        'NumberOfReadsUnmappedTooShort': ["1166"],
        'PectOfReadsUnmappedTooShort': ["87.87"],
        'NumberOfReadsUnmappedOther': ["2"],
        'PectOfReadsUnmappedOther': ["0.15"],
        'ReadsWithValidBarcodes': ["0.023361"],
        'SequencingSaturation': [pd.NA],
        'Q30BasesInCB+UMI': ["0.917408"],
        'ReadsMappedToTranscriptome:Unique+MultipeGenes': ["0"],
        'EstimatedNumberOfCells': ["0"],
        'FractionOfReadsInCells': [pd.NA],
        'MeanReadsPerCell': ["0"],
        'NumberOfUMIs': ["0"],
        'NumberOfGenes': ["0"],
        'NumberOfCountedReads': ["0"],
    }
    expected = pd.DataFrame.from_dict(expected_dict, dtype=pd.StringDtype())
    expected.index = pd.Index(["foo"], name="WellBC", dtype=pd.StringDtype())
    contents = pd.read_csv(output_path, sep="\t", index_col=0, dtype=pd.StringDtype())
    pd.testing.assert_frame_equal(contents, expected)



if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))