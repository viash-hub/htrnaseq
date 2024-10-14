import pytest
import sys
import pandas as pd
from pathlib import Path
from uuid import uuid4

### VIASH START
meta = {
    "resources_dir": "./src/eset/create_pdata/",
    "executable": "target/executable/eset/create_pdata/create_pdata",
    "config": "src/eset/create_pdata/config.vsh.yaml"
}
### VIASH END

@pytest.fixture
def test_reads_and_genes_per_chr_path():
    return Path(meta["resources_dir"]) / "nrReadsNrGenesPerChromPool.txt"


@pytest.fixture
def test_star_logs_summary_path():
    return Path(meta["resources_dir"]) / "starLogs.txt"


@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper 


def test_create_fdata(run_component, test_reads_and_genes_per_chr_path,
                      test_star_logs_summary_path, random_path):
    output_path = random_path("tsv")
    run_component([
        "--star_stats_file", test_star_logs_summary_path,
        "--nrReadsNrGenesPerChromPool", test_reads_and_genes_per_chr_path, 
        "--output", output_path
    ])
    assert output_path.is_file()
    result = pd.read_csv(output_path, sep="\t", dtype=pd.StringDtype())
    expected_dict = {
        'WellBC': ['AACAAGGTAC', 'ACGCCTTCGT', 'CCATACTGAC', 'GCAAGCGAAT',
                   'GTCTCGAGTG', 'TGCGCTCATT', 'TTGTGTTCGA'],
        'NumberOfMTReads': ['0', '0', '0', '0', '0', '0', '0'],
        'pctMT': ['0', '0', '0', '0', '0', '0', '0'],
        'NumberOfERCCReads': ['0', '0', '0', '0', '0', '0', '0'],
        'pctERCC': ['0', '0', '0', '0', '0', '0', '0'],
        'NumberOfChromReads': ['8542', '5863', '7396', '10092', '470',
                               '7650', '9422'],
        'pctChrom': ['100', '100', '100', '100', '100', '100', '100'],
        'NumberOfInputReads': ['141303', '96430', '113577', '156134', '10158',
                               '126989', '142560'],
        'NumberOfMappedReads': ['23749', '16869', '17319', '24005', '1902',
                                '19272', '22129'],
        'PctMappedReads': ['16.81', '17.49', '15.25', '15.37', '18.72',
                           '15.18', '15.52'],
        'NumberOfReadsMappedToMultipleLoci': ['0', '0', '0', '0', '0', '0', '0'],
        'PectOfReadsMappedToMultipleLoci': ['0', '0', '0', '0', '0', '0', '0'],
        'NumberOfReadsMappedToTooManyLoci': ['8458', '6124', '5905', '7961', '967',
                                             '7141', '7045'],
        'PectOfReadsMappedToTooManyLoci': ['5.99', '6.35', '5.2', '5.1', '9.52',
                                           '5.62', '4.94'],
        'NumberOfReadsUnmappedTooManyMismatches': ['0', '0', '0', '0', '0', '0', '0'],
        'PectOfReadsUnmappedTooManyMismatches': ['0', '0', '0', '0', '0', '0', '0'],
        'NumberOfReadsUnmappedTooShort': ['109035', '73375', '90292', '124096',
                                          '7280', '100515', '113324'],
        'PectOfReadsUnmappedTooShort': ['77.16', '76.09', '79.5', '79.48',
                                        '71.67', '79.15', '79.49'],
        'NumberOfReadsUnmappedOther': ['61', '62', '61', '72', '9', '61', '62'],
        'PectOfReadsUnmappedOther': ['0.04', '0.06', '0.05', '0.05',
                                     '0.09', '0.05', '0.04'],
        'ReadsWithValidBarcodes': ['0.999816', '0.999782', '0.999859', '0.999744',
                                   '0.999803', '0.999843', '0.999783'],
        'SequencingSaturation': ['0.0698056', '0.0665302', '0.0717282', '0.0680872',
                                 '0.0553191', '0.0667974', '0.060828'],
        'Q30BasesInCB+UMI': ['0.979965', '0.980077', '0.982313', '0.982779',
                             '0.984451', '0.986581', '0.986622'],
        'ReadsMappedToTranscriptome:Unique+MultipeGenes': ['0.0618175', '0.0620969',
                                                           '0.066554', '0.0658665',
                                                           '0.0476472', '0.0616668',
                                                           '0.0676838'],
        'EstimatedNumberOfCells': ['1', '1', '1', '1', '1', '1', '1'],
        'FractionOfReadsInCells': ['1', '1', '1', '1', '1', '1', '1'],
        'MeanReadsPerCell': ['8538', '5862', '7389',
                             '10090', '470', '7650', '9420'],
        'NumberOfUMIs': ['7942', '5472', '6859', '9403',
                         '444', '7139', '8847'],
        'NumberOfGenes': ['408', '377', '391', '420', '150', '407', '420'],
        'NumberOfCountedReads': ['9535', '6463', '8299', '11273',
                                 '533', '8444', '10383']
    }
    expected = pd.DataFrame.from_dict(expected_dict, dtype=pd.StringDtype())
    pd.testing.assert_frame_equal(result, expected, check_like=True)

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))