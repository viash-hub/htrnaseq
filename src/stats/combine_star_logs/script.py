import logging
import pandas as pd
import numpy as np
from itertools import batched, starmap

### VIASH START
meta = {
    "name": "combine_star_logs",
}
par = {
    "star_logs": ["testData/STAR/ACGCCTTCGT/Log.final.out",
                  "testData/STAR/GTCTCGAGTG/Log.final.out"],
    "gene_summary_logs": ["testData/STAR/ACGCCTTCGT/Solo.out/Gene/Summary.csv",
                          "testData/STAR/GTCTCGAGTG/Solo.out/Gene/Summary.csv"], 
    "reads_per_gene_logs": ["testData/STAR/ACGCCTTCGT/ReadsPerGene.out.tab",
                            "testData/STAR/GTCTCGAGTG/ReadsPerGene.out.tab"],
    "output": "output.txt",
    "barcodes": ["ACGG", "TTTT"],
}

### VIASH END

logger = logging.getLogger()
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)


def handle_percentages(column_value):
    # TODO: handle this more gracefully
    if column_value:
        return np.float64(column_value.strip('%'))
    return column_value

def star_log_to_dataframe(barcode: str, log_path) -> pd.DataFrame:
    logger.info("Reading STAR log %s for barcode '%s'", log_path, barcode)
    result = pd.read_table(log_path, sep=r"\|\t+", converters={"Value": handle_percentages},
                           engine="python", header=None, skip_blank_lines=True,
                           skipinitialspace=True, names=["Category", "Value"], index_col=0,
                           skiprows=[0, 1, 2])
    logger.info("Read %d row(s) and %d column(s) from STAR logs at %s", 
                *result.shape, log_path)
    return result


def summary_to_dataframe(barcode: str, summary_path) -> pd.DataFrame:
    logger.info("Reading summary log %s for barcode %s", summary_path, barcode)
    result = pd.read_table(summary_path, sep=",",
                           header=None, names=["Category", "Value"],
                           index_col=0)
    logger.info("Read %d row(s) and %d column(s) from summary file at %s",
                *result.shape, summary_path)
    return result


def reads_per_gene_to_dataframe(barcode, read_per_gene_path) -> pd.DataFrame:
    logger.info("Reading reads per gene file %s for barcode %s", read_per_gene_path, barcode)
    result = pd.read_table(read_per_gene_path, skiprows=[0, 1, 2, 3], header=None, sep="\t",
                           index_col=0, names=["geneID", "Unstranded", "posStrand", "negStrand"])
    result = result[["Unstranded"]] # Do not use .loc here because we need a DataFrame, not a Series
    df = pd.DataFrame({"Value": result.sum()})
    df = df.rename({"Unstranded": "NumberOfCountedReads"}, errors="raise")
    df.index.name = "Category"
    logger.info("Read %d row(s) and %d column(s) from reads per gene file at %s",
                *df.shape, read_per_gene_path)
    return df

def star_log_remove_unwanted_entries_and_adjust_format(barcode, df: pd.DataFrame) -> pd.DataFrame:
    """
    For a single star log (Log.final.out) in dataframe format, filter out the
    entries that are not needed and format the labels for some metrics:
        - Replace '%' with 'pect' in the labels.
        - Remove labels ending with ':' 
          (mostly the section separators like 'MULTI-MAPPING READS:' and 'UNMAPPED READS:')
        - Remove the metrics we do no need based on the following keywords:
          Mapping speed, Average, Number of splices, per base, chimeric reads, average
    
    The dataframe provided as input must have an index with 1 level with the metric names.
    """
    # Remove index values ending with ':' (rows like 'MULTI-MAPPING READS:','UNIQUE READS:')
    logger.info("Filtering STAR logs for barcode %s. Starting with %d row(s) and %d column(s)", barcode, *df.shape)
    to_keep = ~df.index.to_series().str.endswith(":")
    # Remove index values where the values contain any of these substrings
    regex_columns_to_remove = "Mapping speed|Average|Number of splices|per base|chimeric reads|average"
    to_keep = to_keep & ~df.index.to_series().str.contains(regex_columns_to_remove, regex=True)
    logger.info("Removed the following log entries for barcode '%s':\n\t%s",
                barcode,
                "\n\t".join(to_keep[~to_keep].index.to_list()))
    result = df.loc[to_keep]

    # Replace % by pect, remove columns, use camel case and remove spaces
    # You might be tempted to use .title() to make everything uppercase,
    # but characters which are already uppercase should stay that way.
    # (example: NumberOfUMIs and not NumberOfUmis)
    result.index = result.index.str.replace("%", "pect")\
                    .str.replace(":", "")\
                    .str.replace(r"(?:^|\s).", lambda m:m.group(0).upper(), regex=True)\
                    .str.replace(" ", "")
    result = result.rename({"UniquelyMappedReadsNumber": "NumberOfMappedReads", 
                            "UniquelyMappedReadsPect": "PctMappedReads"}, errors="raise")
    logger.info("Done filtering STAR logs for barcode %s. Result has %d row(s) and %d column(s). "
                "Found entries:\n\t%s", 
                barcode, *result.shape, "\n\t".join(result.index.to_list()))
    return result


def summary_remove_unwanted_entries_and_adjust_format(barcode, df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Filtering and formatting summary logs for barcode %s. "
                "Starting with %d row(s) and %d column(s)", barcode, *df.shape)
    columns_to_remove = (
        "Number of Reads",
        "Q30 Bases in RNA read",
        "Reads Mapped to Genome: Unique",
        "Reads Mapped to Transcriptome: Unique Genes",
        "Reads in Cells Mapped to Unique Genes",
        "Mean Reads per Cell",
        "Median UMI per Cell",
        "Median Genes per Cell",
        "Q30 Bases in CB+UMI",
        "Reads Mapped to Genome: Unique+Multiple",
        "Reads Mapped to Transcriptome: Unique+Multipe Genes",
        "Fraction of Reads in Cells",
        "Median Reads per Cell",
        "Mean UMI per Cell",
        "Mean Genes per Cell",
    )

    to_keep = ~df.index.isin(columns_to_remove)
    logger.info("Removed the following summary entries for barcode '%s':\n\t%s",
                barcode,
                "\n\t".join(df.loc[~to_keep].index.to_list()))
    result = df.loc[to_keep]
    result.index = result.index.str.replace(r"(?:^|\s).", lambda m:m.group(0).upper(),
                                            regex=True).str.replace(" ", "")
    to_rename = {"UMIsInCells": "NumberOfUMIs", 
                 "TotalGenesDetected": "NumberOfGenes"}
    try:
        result = result.rename(to_rename, errors="raise")
    except KeyError as e:
        raise KeyError(f"Tried to rename log entries ({','.join(to_rename)}) in the summary "
                       f"log for barcode {barcode}, but an entry was not found in the file. "
                       "Make sure that you are using the correct version of STAR."
                       f"Available entries: {", ".join(result.index.to_list())}") from e
    logger.info("Done filtering summary logs for barcode %s. Result has %d row(s) and %d column(s). "
                "Found entries:\n\t%s",
                barcode, *result.shape, "\n\t".join(result.index.to_list()))
    return result


def join_dfs(df_list, barcodes) -> pd.DataFrame:
    # Combine the dataframes together and add the barcodes as a level to the dataframe
    # in order to make a 2-level index (first level the barcodes and second level the metrics).
    result = pd.concat(dict(zip(barcodes, df_list)), names=["WellBC"])
    # Pivot the table by moving the metrics to the columns. Its added as an extra level, 
    # so we can just frop the 'Values' level that was already there
    result = result.unstack(level="Category").droplevel(0, axis="columns")
    return result

def main(par):
    logger.info("Component started.")
    # Provide an overview of the parameters in the logs
    parameters_str = [f'\t{param}: {param_val}\n' for param, param_val in par.items()]
    logger.info("Parameters:\n%s", "".join(parameters_str).rstrip())
    star_logs, gene_summary_logs, reads_per_gene_logs, barcodes  = par["star_logs"], \
        par["gene_summary_logs"], par["reads_per_gene_logs"], par["barcodes"]
    number_of_inputs = tuple(len(i) for i in (star_logs, gene_summary_logs,
                                              reads_per_gene_logs, barcodes))
    if len(set(number_of_inputs)) != 1:
        raise ValueError("Expected the same number of inputs for 'star_logs' (%d), "
                         "'gene_summary_logs' (%d), 'reads_per_gene_logs' (%d) "
                         "and 'barcodes' (%d)." % number_of_inputs)
    
    logs_to_process = [
        (star_log_to_dataframe, star_log_remove_unwanted_entries_and_adjust_format, star_logs),
        (summary_to_dataframe, summary_remove_unwanted_entries_and_adjust_format, gene_summary_logs),
        (reads_per_gene_to_dataframe, None, reads_per_gene_logs),
    ]
    logger.info("Formatting the contents of the log files.") 
    all_logs_data = []
    for df_generator, formatter, data in logs_to_process:
        data_as_df = list(starmap(df_generator, zip(barcodes, data)))
        data_formatted = data_as_df
        if formatter:
            data_formatted = list(starmap(formatter, zip(barcodes, data_as_df)))
        data_joined = join_dfs(data_formatted, barcodes)
        all_logs_data.append(data_joined)

    logger.info("Joining entries across the different logs together.") 
    all_stats = pd.concat(all_logs_data, axis=1)
    logger.info("Log statistics were gathered for the following barcodes: %s", 
                ", ".join(all_stats.index.to_list()))
    # batched() is used here to print a limited amount of columnns at a time
    # to make sure that they are all displayed (pandas might limit the view for readability)
    logger.info("Summary of final output:\n%s\n",
                "\n".join(repr(all_stats.loc[:,columns].describe())
                          for columns in batched(all_stats.columns, 3))) 
    logger.info("Writing output to %s", par["output"])
    all_stats.reset_index("WellBC").to_csv(par["output"], sep="\t", header=True, index=False)
    logger.info("Finished %s.", meta["name"])

if __name__ == "__main__":
    main(par)