import pysam
import pandas as pd

## VIASH STAR
par = {


}
## VIASH END

if __name__ == "__main__":
    samfile = pysam.AlignmentFile(par["input"], "rb", threads=2)
    all_tags = []
    index = []
    tags_selection = ("CB", "UB", "GX", "GN")
    for aligned_segment in samfile:
        tags = dict(aligned_segment.get_tags())
        all_tags.append(tags)
        reference_name = aligned_segment.reference_name
        reference_name = "*" if not reference_name else reference_name
        index.append(reference_name)
    tag_dataframe = pd.DataFrame.from_records(all_tags, index=index)
    
    tag_dataframe.to_csv(par["processedBamFile"], sep="\t", na_rep="",
                         header=True, index_label="Chr",
                         columns=["CB", "UB", "GX", "GN"])

    # Number of genes that had a read mapped to them per chromosome,
    # and the number of reads mapped to those genes pet chromosome.
    nr_reads_nr_genes = tag_dataframe.dropna(subset=["GX"]).groupby(level=0).agg(
        NumberOfReads=pd.NamedAgg("GX", aggfunc="size"),
        NumberOfGenes=pd.NamedAgg(column="GX", aggfunc="nunique")
    )
    nr_reads_nr_genes.to_csv(par["nrReadsNrGenesPerChrom"], sep="\t", header=True,
                             index_label="Chr") 
    
    # Number of reads mapped to the reference, grouped by UMI
    nr_read_per_umi = tag_dataframe.groupby('UB').size()\
        .drop("", errors="ignore").sort_values(ascending=False).head(100) 
    nr_read_per_umi.to_frame(name="N").to_csv(par["umiFreqTop100"], header=True, sep="\t")

    # Total number of mapped reads and total number of UMIs (not grouped per chromosome)
    nr_reads_and_umi_per_barcode = tag_dataframe.groupby(by="CB").agg(
        NumberOfReads=pd.NamedAgg("CB", "size"),
        nrUMIs=pd.NamedAgg("UB", "nunique")
    )
    nr_reads_and_umi_per_barcode.to_csv(par["nrReadsNrUMIsPerCB"], sep="\t", header=True)