from itertools import batched
import pandas as pd
import logging

### VIASH START
meta = {
    "name": "create_pdata",
}

par = {
  "star_stats_file": "src/eset/create_pdata/starLogs.txt",
  "nrReadsNrGenesPerChromPool": "src/eset/create_pdata/nrReadsNrGenesPerChromPool.txt",
  "output": "pData.tsv"
}

### VIASH END

logger = logging.getLogger()
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

def main(par):
  logger.info(f"{meta['name']} started.")
  parameters_str = [f'\t{param}: {param_val}\n' for param, param_val in par.items()]
  logger.info("Parameters:\n%s", "".join(parameters_str).rstrip())
  logger.info("Reading %s", par["star_stats_file"])
  star_log_stats = pd.read_csv(par["star_stats_file"], sep="\t", index_col=0)
  logger.info("STAR log statics file contains information for the following barcodes: %s", 
              ", ".join(star_log_stats.index))
  logger.info("Reading %s", par["nrReadsNrGenesPerChromPool"])
  reads_and_genes_per_chr_stats = pd.read_csv(par["nrReadsNrGenesPerChromPool"], sep="\t", index_col=0)
  logger.info("Reads per gene and chromosome table contains information for the following barcodes: %s",
              ", ".join(reads_and_genes_per_chr_stats.index))
  logger.info("Filtering mapping statistics file columns.")
  cols_to_keep = ("WellID", "NumberOfMTReads", "pctMT", "NumberOfERCCReads",
                  "pctERCC", "NumberOfChromReads", "pctChrom")
  try:
    reads_and_genes_per_chr_stats = reads_and_genes_per_chr_stats.loc[:,cols_to_keep]
  except KeyError as e:
    raise KeyError("When trying to subset the reads per genes and chromosomes file, "
                   "a column was missing. Available columns in the file: "
                   f"{', '.join(reads_and_genes_per_chr_stats.columns)}.") from e
  # Each barcode should be present. An alternative approach could be to just
  # do the concatenation and check for NA values that are filled for non-overlapping
  # index values, but there are already NA values present in the dataframes
  if not star_log_stats.index.sort_values().equals(reads_and_genes_per_chr_stats.index.sort_values()):
    raise ValueError("Error while combining two log files. It seems that the entries (barcodes) "
                     f"do not fully overlap. Barcodes in '{par['star_stats_file']}: "
                     f"{', '.join(reads_and_genes_per_chr_stats.index)}. Barcodes in "
                     f"'{par['nrReadsNrGenesPerChromPool']}': "
                     f"{', '.join(star_log_stats.index)}")
  combined_stats = pd.concat([reads_and_genes_per_chr_stats, star_log_stats], axis=1)
  logger.info("Summary of final output:\n%s\n",
                "\n".join(repr(combined_stats.loc[:,columns].describe())
                          for columns in batched(combined_stats.columns, 3))) 
  logger.info("Writing to %s", par["output"])
  combined_stats.reset_index("WellBC").to_csv(par["output"], sep="\t", header=True, index=False)
  logger.info("Finished %s.", meta["name"])


if __name__ == "__main__":
  main(par)