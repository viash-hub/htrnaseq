import pysam
import pandas as pd
import logging

### VIASH START
par = {
    "input": "src/stats/generate_well_statistics/test.sam",
    "processedBAMFile": "processedBamFile.txt",
    "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom.txt",
    "nrReadsNrUMIsPerCB": "nrReadsNrUMIsPerCB.txt",
    "umiFreqTop": "umiFreqTop.txt",
    "threads": 1,
    "barcode": "ACGT"
}
### VIASH END
logger = logging.getLogger()
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

if __name__ == "__main__":
    logger.info("Component started.")
    parameters_str = [f'\t{param}: {param_val}\n' for param, param_val in par.items()]
    logger.info("Parameters:\n%s", "".join(parameters_str).rstrip())
    logger.info("Opening '%s'", par["input"])
    samfile = pysam.AlignmentFile(par["input"], "rb", threads=par["threads"])
    all_tags = []
    index = []
    tags_selection = ("CB", "UB", "GX", "GN")
    for aligned_segment in samfile:
        tags = dict(aligned_segment.get_tags())
        all_tags.append(tags)
        reference_name = aligned_segment.reference_name
        index.append("*" if not reference_name else reference_name)
    tag_dataframe = pd.DataFrame.from_records(all_tags, index=index,
                                              columns=tags_selection)
    tag_dataframe_to_write = tag_dataframe.copy()
    logger.info("Done reading BAM file. Found %i entries", tag_dataframe.shape[0])
    tag_dataframe.assign(WellBC=par["barcode"], WellID=par["well_id"])\
        .reset_index(names="Chr")\
        .to_csv(par["processedBAMFile"], sep="\t", na_rep="",
                header=True, index=False,
                columns=("WellBC", "WellID", "Chr") + tags_selection)
    logger.info("Constructing of dataframe done.")
    # Number of genes that had a read mapped to them per chromosome,
    # and the number of reads mapped to those genes per chromosome.
    nr_reads_nr_genes = tag_dataframe.dropna(subset=["GX"]).groupby(level=0).agg(
        NumberOfReads=pd.NamedAgg("GX", aggfunc="size"),
        NumberOfGenes=pd.NamedAgg(column="GX", aggfunc="nunique")
    )
    logger.info("Done calculating number of reads per gene and per chromesome. Writing to %s",
                par['nrReadsNrGenesPerChrom'])
    nr_reads_nr_genes.reset_index(names="Chr").assign(WellBC=par["barcode"], WellID=par["well_id"])\
        .to_csv(par["nrReadsNrGenesPerChrom"], sep="\t",
                header=True, index=False, 
                columns=("WellBC", "WellID", "Chr", "NumberOfReads", "NumberOfGenes"))

    # Number of reads mapped to the reference, grouped by UMI
    nr_read_per_umi = tag_dataframe.groupby('UB').size()\
        .drop("", errors="ignore").sort_values(ascending=False).head(100)
    nr_read_per_umi_df = nr_read_per_umi.to_frame(name="N")
    logger.info("Done calculating number of mapped reads per UMI, writing to %s", par["umiFreqTop"])
    nr_read_per_umi_df.assign(WellBC=par["barcode"], WellID=par["well_id"]).reset_index(names="UB")\
        .to_csv(par["umiFreqTop"], header=True, sep="\t", 
                index=False, columns=("WellBC", "WellID", "UB", "N"))

    # Total number of mapped reads and total number of UMIs (not grouped per chromosome)
    nr_reads_and_umi_per_barcode = tag_dataframe.groupby(by="CB").agg(
        NumberOfReads=pd.NamedAgg("CB", "size"),
        nrUMIs=pd.NamedAgg("UB", "nunique")
    )
    logger.info("Done calculating number of mapped reads and number of UMIs per Cell Barcode, writing to %s",
                par["nrReadsNrUMIsPerCB"])
    nr_reads_and_umi_per_barcode.assign(WellBC=par["barcode"], WellID=par["well_id"]).reset_index(names="CB")\
        .to_csv(par["nrReadsNrUMIsPerCB"], sep="\t", header=True, 
                index=False, columns=("WellBC", "WellID", "CB", "NumberOfReads", "nrUMIs"))
    logger.info("Finished!")