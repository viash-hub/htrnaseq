import logging
import pandas as pd
import numpy as np
from textwrap import fill


### VIASH START
meta = {
    "name": "create_fdata",
}

par = {
  "gtf": "src/eset/create_fdata/Homo_sapiens_ERCC.rnaseq_starReference_0.0.3.chr.gtf",
  "output": "fData.tsv"
}

### VIASH END

logger = logging.getLogger()
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)


def read_gtf(gtf_path: str) -> pd.DataFrame:
    logger.info("Reading %s", gtf_path)
    result = pd.read_csv(gtf_path, sep="\t",
                         header=None, names=("seqname", "source",
                                             "feature", "start", "end",
                                             "score", "strand", "frame",
                                             "attribute"),
                         dtype={
                            "seqname": pd.StringDtype(),
                            "source": pd.StringDtype(),
                            "feature": pd.StringDtype(),
                            "start": pd.Int64Dtype(),
                            "end": pd.Int64Dtype(),
                            "score": pd.StringDtype(),
                            "strand": pd.CategoricalDtype(categories=["+", "-"],
                                                            ordered=False),
                            "frame": pd.StringDtype(),
                            "attribute": pd.StringDtype(),
                          },
                          comment='#'
                        )
    logger.info("Done reading %s. Found %d GTF entries ", par["gtf"], result.shape[0])
    logger.info("GTF file is providing information for the following chromosomes: \n%s", 
                fill(", ".join(result['seqname'].unique()), width=100))
    logger.info("The following sources were specified in the GTF file:\n%s",
                ", ".join(result["source"].unique()))
    return result
    

def parse_attributes(attributes_series: pd.Series):
    attribute_dict = dict()
    attributes_list = [attr.strip().split(" ")
                       for attr in attributes_series["attribute"].strip(";").split(";")]
    for (attr_name, attr_value) in attributes_list:
        attribute_dict.setdefault(attr_name, []).append(attr_value.strip('"'))
    attribute_dict = {attr_name: "|".join(attr_value) 
                      for attr_name, attr_value in attribute_dict.items()}
    return pd.Series(attribute_dict)
    

def main(par):
    logger.info(f"{meta['name']} started.")
    parameters_str = [f'\t{param}: {param_val}\n' for param, param_val in par.items()]
    logger.info("Parameters:\n%s", "".join(parameters_str).rstrip())
    gtf_file = read_gtf(par["gtf"])
    sources = set(source for source in gtf_file["source"].unique() if source != "ERCC")
    specific_gtf = False
    feature = "gene"
    if len(sources) == 1 and (source := sources[0]) \
        and (source == "refGene" or source == "ncbiRefSeq"):
        feature = "transcript"
        specific_gtf = True
        logger.info("Found specific GTF from %s, forcing filtering on feature type %s", source, feature)
    logger.info("Filtering GTF entries for feature type '%s'.", feature)
    gtf_file = gtf_file[gtf_file["feature"] == feature]
    logger.info("After filtering %d entries are left.", gtf_file.shape[0])
    logger.info("Parsing the GTF attributes")
    annotation = gtf_file[["attribute"]].apply(parse_attributes, result_type="expand", axis=1)
    logger.info("Found the following attributes in the GTF:\n%s", ", ".join(annotation.columns))
    annotation = pd.concat([gtf_file.drop(["attribute"], axis=1), annotation], axis=1)
    if specific_gtf:
       logger.info("Because the source of the GTF is either 'ncbiRefSeq' or 'refGene', which"
                   "caused forced filtering based on %s, the duplicate genes still need to be dropped.",
                   feature)
       annotation = annotation.drop_duplicates(subset=("gene_id", "gene_name"), keep=False)
       logger.info("After dropping duplicates, %d entries are left", annotation.shape[0])

    # detect ensembl ids
    # some GTF files contain version in ENSEMBL, e.g. ENS00000000046319.1
    # we remove the version, because the annotation packages don't contain the version
    if "gene_id" in annotation.columns:
        logger.info("'gene_id' column was detected in attributes. Performing extra parsing of ENSEMBL ids.")
        annotation["ENSEMBL_with_version"] = annotation["gene_id"].where(annotation["gene_id"].str.startswith("ENS"))
        annotation["ENSEMBL"] = annotation["ENSEMBL_with_version"].str.replace(r"\.\d+$", "", regex=True)
        annotation["gene_id"] = annotation["gene_id"].str.replace(r"\.\d+$", "", regex=True)

    possible_name_columns = ("Name", "name", "gene_name")
    found_columns = list(filter(lambda col_name: col_name in annotation, possible_name_columns))
    # The following code allows to select a value for the SYMBOL column based on the first non-na column
    if found_columns:
        logger.info("Found one the following columns: %s; which can be used to populate the SYMBOL column",
                    ", ".join(possible_name_columns))
        # For each row (gtf entry), get the name of the first column that actually holds a value.
        column_to_get = annotation.loc[:,found_columns].apply(pd.Series.first_valid_index, axis=1)
        counts_per_column = column_to_get.value_counts(dropna=False).to_dict()
        counts_per_column_str = [f'\t{col}: {counts}\n' for col, counts in counts_per_column.items()]
        logger.info("Frequencies of the origin for the entries in the SYMBOL column:\n%s",
                    "".join(counts_per_column_str).rstrip())
        # If all columns hold NA for a certain row, first_valid_index will return None.
        # Just use the name of the first column.
        column_to_get = column_to_get.fillna(found_columns[0])
        # We now have a list one column name per row, use it so select the values
        # Loc cannot be used here because 1 value per row is required, 
        # and loc will select for each row all the columns in columns_to_get
        idx, cols = pd.factorize(column_to_get)
        symbol_values = annotation.reindex(cols, axis=1).to_numpy()[np.arange(len(annotation)), idx]
        annotation["SYMBOL"] = symbol_values
    logger.info("Dropping unused columns")
    annotation = annotation.drop(["score", "source", "frame", "feature"], axis=1)
    logger.info("Looking for duplicate rows and removing them. Starting with %i entries", annotation.shape[0])
    annotation = annotation.drop_duplicates(keep="first", ignore_index=True)
    logger.info("After removing duplicates: %i entries", annotation.shape[0])
    logger.info("Writing to %s", par["output"])
    annotation.to_csv(par["output"], sep="\t", header=True, index=False, na_rep="NA")
    # Do these checks *after* writing the csv in order to be able to check the data
    logger.info("Checking for unique gene IDs")
    if not annotation["gene_id"].is_unique:
        raise ValueError("Values from the 'gene_id' column are not unique after processing!") 
    logger.info("%s finished", meta['name'])


if __name__ == "__main__":
    main(par)
