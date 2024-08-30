import logging
import pandas as pd
import numpy as np

### VIASH START
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

def star_log_to_dataframe(log_path):
    return pd.read_table(log_path, sep=r"\|\t+", converters={"Value": handle_percentages},
                         engine="python", header=None, skip_blank_lines=True,
                         skipinitialspace=True, names=["Category", "Value"], index_col=0,
                         skiprows=[0, 1, 2])

def summary_to_dataframe(summary_path):
    return pd.read_table(summary_path, sep=",",
                         header=None, names=["Category", "Value"],
                         index_col=0)


def reads_per_gene_to_dataframe(read_per_gene_path):
    result = pd.read_table(read_per_gene_path, skiprows=[0, 1, 2, 3], header=None, sep="\t",
                           index_col=0, names=["geneID", "Unstranded", "posStrand", "negStrand"])
    result = result[["Unstranded"]] # Do not use .loc here because we need a DataFrame, not a Series
    df = pd.DataFrame({"Value": result.sum()})
    df.index.name = "Category"
    return df

def star_log_remove_unwanted_entries_and_adjust_format(df: pd.DataFrame) -> pd.DataFrame:
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
    to_keep = ~df.index.to_series().str.endswith(":")
    # Remove index values where the values contain any of these substrings
    regex_columns_to_remove = "Mapping speed|Average|Number of splices|per base|chimeric reads|average"
    to_keep = to_keep & ~df.index.to_series().str.contains(regex_columns_to_remove, regex=True)
    result = df.loc[to_keep]

    result.index = result.index.str.replace("%", "pect")\
                    .str.replace(":", "")\
                    .str.replace(r"(?:^|\s).", lambda m:m.group(0).upper(), regex=True)\
                    .str.replace(" ", "")
    result = result.rename({"UniquelyMappedReadsNumber": "NumberOfMappedReads", 
                            "UniquelyMappedReadsPect": "pctMappedReads"}, errors="raise")
    return result


def summary_remove_unwanted_entries_and_adjust_format(df: pd.DataFrame) -> pd.DataFrame:
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
    result = df.loc[to_keep]
    result.index = result.index.str.replace(r"(?:^|\s).", lambda m:m.group(0).upper(), regex=True).str.replace(" ", "")
    result = result.rename({"UMIsInCells": "NumberOfUMIs", 
                            "TotalGenesDetected": "NumberOfGenes"}, errors="raise")
    return result


def join_dfs(df_list, barcodes):
    # Combine the dataframes together and add the barcodes as a level to the dataframe
    # in order to make a 2-level index (first level the barcodes and second level the metrics).
    result = pd.concat(dict(zip(barcodes, df_list)), names=["WellBC"])
    # Pivot the table by moving the metrics to the columns. Its added as an extra level, 
    # so we can just frop the 'Values' level that was already there
    result = result.unstack(level="Category").droplevel(0, axis="columns")
    return result

def main(par):
    logger.info("Component started.")
    parameters_str = [f'\t{param}: {param_val}\n' for param, param_val in par.items()]
    logger.info("Parameters:\n%s", "".join(parameters_str).rstrip())
    star_logs, gene_summary_logs, reads_per_gene_logs, barcodes  = par["star_logs"], \
        par["gene_summary_logs"], par["reads_per_gene_logs"], par["barcodes"]
    
    if len({len(i) for i in (star_logs, gene_summary_logs, reads_per_gene_logs,
                             reads_per_gene_logs, barcodes)}) != 1:
        raise ValueError("Expected the same number of inputs for 'star_logs', "
                         "'gene_summary_logs' and 'reads_per_gene_logs'")
    
    # TODO: put this in a loop
    star_log_dfs = list(map(star_log_to_dataframe, star_logs))
    star_log_dfs = list(map(star_log_remove_unwanted_entries_and_adjust_format, star_log_dfs))
    star_log_df = join_dfs(star_log_dfs, barcodes)

    summary_dfs = list(map(summary_to_dataframe, gene_summary_logs))
    summary_dfs = list(map(summary_remove_unwanted_entries_and_adjust_format, summary_dfs))
    summary_df = join_dfs(summary_dfs, barcodes)

    gene_tab_dfs = list(map(reads_per_gene_to_dataframe, reads_per_gene_logs))
    gene_tab_df = join_dfs(gene_tab_dfs, barcodes)

    all_stats = pd.concat([star_log_df, summary_df, gene_tab_df], axis=1)
    all_stats.reset_index("WellBC").to_csv(par["output"], sep="\t", header=True, index=False)


if __name__ == "__main__":
    main(par)