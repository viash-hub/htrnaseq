import pandas as pd
from pathlib import Path
import re

### VIASH START
par = {
    "nrReadsNrGenesPerChrom": ["src/stats/generate_pool_statistics/test1.tsv", "src/stats/generate_pool_statistics/test2.tsv"],
    "nrReadsNrGenesPerChromPool": "nrReadsNrGenesPerChrom_pool.txt"
}

### VIASH END

INDEX_COL = ["WellBC", "WellID"]

if __name__ == "__main__":
    #########
    # nrReadsNrGenesPerChrom file
    #########
    nr_reads_nr_genes_wells = []
    par["nrReadsNrGenesPerChrom"] = list(map(Path, par["nrReadsNrGenesPerChrom"]))
    for nr_reads_nr_genes_file in par["nrReadsNrGenesPerChrom"]:
        nr_reads_nr_gene_well = pd.read_csv(nr_reads_nr_genes_file,
                                            header=0, delimiter="\t",
                                            dtype={"WellBC": pd.StringDtype(),
                                                   "WellID": pd.StringDtype(),
                                                   "Chr": pd.StringDtype(),
                                                   "NumberOfReads": pd.UInt64Dtype(),
                                                   "NumberOfGenes": pd.UInt64Dtype()})
        if nr_reads_nr_gene_well.empty:
            raise ValueError(f"{nr_reads_nr_genes_file.name} does not seem to contain any information!")
        nr_reads_nr_genes_wells.append(nr_reads_nr_gene_well)
    nr_reads_nr_genes_pool = pd.concat(nr_reads_nr_genes_wells, ignore_index=True,)
    total_nr_reads_per_chromosome = nr_reads_nr_genes_pool.pivot_table(index=INDEX_COL, columns="Chr",
                                                                       values=["NumberOfReads"], fill_value=0,
                                                                       aggfunc="sum").droplevel(0, axis=1)
    total_nr_reads_per_chromosome.columns.name = None
    # Remove scaffolds/chromosomes with no counts
    total_nr_reads_per_chromosome = total_nr_reads_per_chromosome.loc[:, (total_nr_reads_per_chromosome != 0).any(axis=0)]
    ##### Total number of genes from all chromosomes
    total_nr_genes = nr_reads_nr_genes_pool.loc[:, INDEX_COL + ['NumberOfGenes']].groupby(["WellBC", "WellID"]).sum()

    ##### Total counts across (irrespective of chromosome)
    total_sum_of_reads = total_nr_reads_per_chromosome.sum(numeric_only=True, axis=1) 

    ##### Logic to split up chromosome per type
    chromosome_names = total_nr_reads_per_chromosome.columns.to_list()
    chr_regex = re.compile(r"^(chr)?\d+")
    matching_chromosomes = [chr_name for chr_name 
                            in chromosome_names
                            if chr_regex.match(chr_name)]
    sex_chromosome_names = ["X", "Y"]
    mitochondrial_chr_name = "MT"
    # This is logic from the original HT pipeline,
    # only when all of the matched chromosomes start with "chr", the mitochonrial, X and Y
    # chromosomes should also start with 'chr'
    if all(chr_name.startswith("chr") for chr_name in matching_chromosomes):
       sex_chromosome_names += ["chrX", "chrY"]
       mitochondrial_chr_name = "chrM"

    ###### Counts for mitochondrial reads
    try:
        mitochondrial_reads = total_nr_reads_per_chromosome.loc[:,mitochondrial_chr_name]
    except KeyError:
       mitochondrial_reads = 0
    percentage_mitochondrial_reads = round(mitochondrial_reads / total_sum_of_reads * 100, 2)

    ###### Counts for ERCC reads
    total_ercc_reads = total_nr_reads_per_chromosome.filter(regex=r"^ERCC").sum(axis=1)
    percentage_ercc_reads = round(total_ercc_reads / total_sum_of_reads * 100, 2)

    ###### Counts for nuclear chromosomes
    total_chromosomal_reads = total_nr_reads_per_chromosome.loc[:,matching_chromosomes].sum(axis=1)
    percentage_chromosomal_reads = round(total_chromosomal_reads / total_sum_of_reads * 100, 2)

    cols_to_add = {
        "pctChrom": percentage_chromosomal_reads,
        "pctMT": percentage_mitochondrial_reads,
        "pctERCC": percentage_ercc_reads,
        "SumReads": total_sum_of_reads,
        "NumberOfGenes": total_nr_genes,
        "NumberOfERCCReads": total_ercc_reads,
        "NumberOfChromReads": total_chromosomal_reads,
        "NumberOfMTReads": mitochondrial_reads,
    }
    total_nr_reads_per_chromosome = total_nr_reads_per_chromosome.assign(
       **cols_to_add
    )

    total_nr_reads_per_chromosome.reset_index(names=INDEX_COL)\
        .to_csv(par["nrReadsNrGenesPerChromPool"], sep="\t",
                header=True, index=False, float_format="%g",
                columns=tuple(INDEX_COL) + tuple(chromosome_names) + tuple(cols_to_add.keys())
               )

