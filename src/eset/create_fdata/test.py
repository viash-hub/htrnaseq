import pytest
import sys
import pandas as pd
from pathlib import Path
from uuid import uuid4
from shutil import copyfile

### VIASH START
meta = {
    "resources_dir": "./src/eset/create_fdata/",
    "executable": "target/executable/eset/create_fdata/create_fdata",
    "config": "src/eset/create_fdata/config.vsh.yaml"
}
### VIASH END

@pytest.fixture
def test_annotation_path():
    return Path(meta["resources_dir"]) / "test_annotation.gtf"


@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
    return wrapper 


def test_create_fdata(run_component, test_annotation_path, random_path):
    output_path = random_path("tsv")
    run_component([
        "--gtf", test_annotation_path,
        "--output", output_path
    ])
    assert output_path.is_file()
    result = pd.read_csv(output_path, sep="\t", dtype=pd.StringDtype())

    expected_dict = {
        "seqname": ["20", "20", "20", "21"],
        "start": ["87250", "142590", "157454", "297570"],
        "end": ["97094", "145751", "159163", "300321"],
        "strand": ["+", "+", "+", "+"],
        "gene_id": ["ENSG00000178591", "ENSG00000125788",
                    "ENSG00000088782", "ENSG00000247315"],
        "gene_version": ["7", "6", "5", "4"],
        "gene_name": ["DEFB125", "DEFB126", "DEFB127", pd.NA],
        "gene_source": ["ensembl_havana", "ensembl_havana",
                        "ensembl_havana", "havana"],
        "gene_biotype": ["protein_coding", "protein_coding",
                         "protein_coding", "protein_coding"],
        "ENSEMBL_with_version": ["ENSG00000178591.7", "ENSG00000125788",
                                 "ENSG00000088782", "ENSG00000247315"],
        "ENSEMBL": ["ENSG00000178591", "ENSG00000125788",
                    "ENSG00000088782", "ENSG00000247315"],
        "SYMBOL": ["DEFB125", "DEFB126", "DEFB127", pd.NA]
    }
    expected = pd.DataFrame.from_dict(expected_dict, dtype=pd.StringDtype())
    pd.testing.assert_frame_equal(expected, result, check_like=True)


def test_make_unique(run_component, test_annotation_path, random_path):
    gtf_with_duplicate_entry_path = random_path("gtf")
    output_path = random_path("tsv")
    entry_to_add = (
        "\n20	ensembl_havana	gene	87250	97094	.	+	.	gene_id " +
        "\"ENSG00000178591.7\"; gene_version \"7\"; gene_name \"DEFB125\"; " +
        "gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\";\n"
    )
    copyfile(test_annotation_path, gtf_with_duplicate_entry_path)
    with gtf_with_duplicate_entry_path.open("a") as open_gtf:
        open_gtf.write(entry_to_add)
    run_component([
        "--gtf", gtf_with_duplicate_entry_path,
        "--output", output_path
    ])
    assert output_path.is_file()
    result = pd.read_csv(output_path, sep="\t", dtype=pd.StringDtype())
    expected_dict = {
        "seqname": ["20", "20", "20", "21"],
        "start": ["87250", "142590", "157454", "297570"],
        "end": ["97094", "145751", "159163", "300321"],
        "strand": ["+", "+", "+", "+"],
        "gene_id": ["ENSG00000178591", "ENSG00000125788",
                    "ENSG00000088782", "ENSG00000247315"],
        "gene_version": ["7", "6", "5", "4"],
        "gene_name": ["DEFB125", "DEFB126", "DEFB127", pd.NA],
        "gene_source": ["ensembl_havana", "ensembl_havana",
                        "ensembl_havana", "havana"],
        "gene_biotype": ["protein_coding", "protein_coding",
                         "protein_coding", "protein_coding"],
        "ENSEMBL_with_version": ["ENSG00000178591.7", "ENSG00000125788",
                                 "ENSG00000088782", "ENSG00000247315"],
        "ENSEMBL": ["ENSG00000178591", "ENSG00000125788",
                    "ENSG00000088782", "ENSG00000247315"],
        "SYMBOL": ["DEFB125", "DEFB126", "DEFB127", pd.NA]
    }
    expected = pd.DataFrame.from_dict(expected_dict, dtype=pd.StringDtype())
    pd.testing.assert_frame_equal(expected, result, check_like=True)


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))