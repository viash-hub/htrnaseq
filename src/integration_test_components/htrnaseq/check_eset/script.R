library(Biobase)
library(testthat)
library(Matrix)

sample_1_result <- readRDS(par$eset)
expected_sample_names <- c(
  "sample_one_AACAAGGTAC", "sample_one_AACAATCAGG", "sample_one_AACACCTAGT",
  "sample_one_AACAGGCAAT", "sample_one_AACATGGAGA", "sample_one_AACATTACCG",
  "sample_one_AACCAGCCAG", "sample_one_AACCAGTTGA", "sample_one_AACCGCGACT",
  "sample_one_AACCGGAAGG", "sample_one_AACCGGCGTA", "sample_one_AACCTAGTCC",
  "sample_one_AACCTCATAG", "sample_one_AACGTAAGCT", "sample_one_AACTCTACAC",
  "sample_one_AACTGTGTCA", "sample_one_AAGACGGATT", "sample_one_AAGATCGGCG",
  "sample_one_AAGATGTCCA", "sample_one_AAGCATATGG", "sample_one_AAGCGATGTT",
  "sample_one_AAGCGTTCAG", "sample_one_AAGCTCACCT", "sample_one_AAGGCATGCG",
  "sample_one_AAGGTCTGGA", "sample_one_AAGTTAGCGC", "sample_one_AAGTTCCTTG",
  "sample_one_AATACCGGTA", "sample_one_AATAGCCACA", "sample_one_AATCACGCGA",
  "sample_one_AATCCATCTG", "sample_one_AATCCGCTCC", "sample_one_AATCCTACCA",
  "sample_one_AATCGTCCGC", "sample_one_AATGAACACG", "sample_one_AATGACCTTC",
  "sample_one_AATGAGAGCA", "sample_one_AATGTCAGTG", "sample_one_AATTAGGCCG",
  "sample_one_AATTGCGATG", "sample_one_ACAACAGTCG", "sample_one_ACAACCATAC",
  "sample_one_ACAACGGAGC", "sample_one_ACAAGCGCGA", "sample_one_ACACAATCTC",
  "sample_one_ACACAGTGAA", "sample_one_ACACCGAATT", "sample_one_ACACGCAGTA",
  "sample_one_ACACGGTCCT", "sample_one_ACACTTGCTG", "sample_one_ACAGTGCCAA",
  "sample_one_ACATGTGTGC", "sample_one_ACCAGGACCA", "sample_one_ACCATAACAC",
  "sample_one_ACCGAACCGT", "sample_one_ACCGAGAGTC", "sample_one_ACCGGTACAG",
  "sample_one_ACCGTACTTC", "sample_one_ACCTCCGACA", "sample_one_ACCTCTCTCC",
  "sample_one_ACCTGTCCGA", "sample_one_ACCTTATGTG", "sample_one_ACGAATGACA",
  "sample_one_ACGCCTCAAC", "sample_one_ACGCCTTCGT", "sample_one_ACGCTGGATA",
  "sample_one_ACGGTCCGTT", "sample_one_ACGTAGGCAC", "sample_one_ACGTGCTGAT",
  "sample_one_ACTCCAAGCC", "sample_one_ACTGGCGCAT", "sample_one_ACTGGCTTCC",
  "sample_one_ACTTAACTGC", "sample_one_ACTTCATCAC", "sample_one_ACTTCGTTGA",
  "sample_one_ACTTCTCCTG", "sample_one_ACTTGAGGAA", "sample_one_ACTTGTAAGG",
  "sample_one_AGAACCACGG", "sample_one_AGAAGCAATC", "sample_one_AGACCGTTAT",
  "sample_one_AGACTAGCAT", "sample_one_AGAGATGCAG", "sample_one_AGAGCTTACA",
  "sample_one_AGAGTGTAAC", "sample_one_AGAGTTCTGC", "sample_one_AGATAGTGCT",
  "sample_one_AGCAATGCGC", "sample_one_AGCATGTCAT", "sample_one_AGCCACTAGC",
  "sample_one_AGCCAGAATA", "sample_one_AGCCAGCTCT", "sample_one_AGCGATAACG",
  "sample_one_AGCGTACAAT", "sample_one_AGCTATTCCA", "sample_one_AGCTCCTCAG",
  "sample_one_AGGAGGCATA", "sample_one_AGGCGTCTGT", "sample_one_AGTAACTCAC",
  "sample_one_AGTAAGCGTT", "sample_one_AGTCTGTACG", "sample_one_AGTGCAATGT",
  "sample_one_ATAAGGTGCA", "sample_one_ATACACGACA", "sample_one_ATAGGCCATT",
  "sample_one_ATATCCGCAT", "sample_one_ATCAGCACTT", "sample_one_ATCAGCGAGG",
  "sample_one_ATCCAATACG", "sample_one_ATCCGCTGTG", "sample_one_ATCCGTCCAT",
  "sample_one_ATCGACGGCT", "sample_one_ATCGCGATTA", "sample_one_ATCGGTAGGC",
  "sample_one_ATCTAAGGAG", "sample_one_ATGACGGTAA", "sample_one_ATGACTCAGT",
  "sample_one_ATGCACCGGA", "sample_one_ATGCGGACTG", "sample_one_ATGCTTCCTA",
  "sample_one_ATGGACCAAC", "sample_one_ATGGTCTTAG", "sample_one_ATGGTGAGCG",
  "sample_one_ATGTGGAAGC", "sample_one_ATTATCGGAC", "sample_one_ATTCGGAACA",
  "sample_one_CAACAATCCA", "sample_one_CAAGAAGCAT", "sample_one_CAAGATGAGG",
  "sample_one_CAAGCCAACG", "sample_one_CAAGTGGATC", "sample_one_CACAGTTCAT",
  "sample_one_CACGAGTCTG", "sample_one_CACGCTCCAA", "sample_one_CACTGAGCAC",
  "sample_one_CAGATCAATG", "sample_one_CAGTGCTCTT", "sample_one_CAGTTAAGCA",
  "sample_one_CATAGCTATC", "sample_one_CATCACCACC", "sample_one_CATGTACGCC",
  "sample_one_CATTACACTG", "sample_one_CATTCGACGA", "sample_one_CCAACTATGG",
  "sample_one_CCAAGGAGTT", "sample_one_CCAATTGTTC", "sample_one_CCACAAGTGC",
  "sample_one_CCAGCTTAGT", "sample_one_CCATAACTTG", "sample_one_CCATACTGAC",
  "sample_one_CCATAGATCA", "sample_one_CCATGTGCTT", "sample_one_CCATTCAGCG",
  "sample_one_CCGAACAAGC", "sample_one_CCGAACCTAA", "sample_one_CCGAAGACCT",
  "sample_one_CCGAATAGTG", "sample_one_CCGACTTCTC", "sample_one_CCGATCCACT",
  "sample_one_CCGATGATAC", "sample_one_CCGCGTTATG", "sample_one_CCGCTAGCTT",
  "sample_one_CCGGAGTATC", "sample_one_CCGGCCAATT", "sample_one_CCGGTCTCTA",
  "sample_one_CCGTACGATG", "sample_one_CCGTCAGAAC", "sample_one_CCTAGACACG",
  "sample_one_CCTAGTTGAG", "sample_one_CCTATTCTGT", "sample_one_CCTCAACCGA",
  "sample_one_CCTCCATAAG", "sample_one_CCTGATGCCA", "sample_one_CCTGCAATAC",
  "sample_one_CCTTGTATTC", "sample_one_CGAGATCTCT", "sample_one_CGAGGAACAA",
  "sample_one_CGATAACCGC", "sample_one_CGATCCTGTG", "sample_one_CGCCAACCAT",
  "sample_one_CGCCAGTGTT", "sample_one_CGCCTTGTAC", "sample_one_CGCGGATTCA",
  "sample_one_CGCTTAAGGC", "sample_one_CGCTTACTAA", "sample_one_CGCTTCTTGG",
  "sample_one_CGGAAGCTGT", "sample_one_CGGAATACAC", "sample_one_CGGAGATTGG",
  "sample_one_CGGAGCTCAA", "sample_one_CGGATCGGTA", "sample_one_CGGATTCTAG",
  "sample_one_CGGCAACTTA", "sample_one_CGGCTCATCA", "sample_one_CGGTCGTATT",
  "sample_one_CGGTGACATC", "sample_one_CGTAACGGAT", "sample_one_CGTAAGATTC",
  "sample_one_CGTACTGTAA", "sample_one_CGTAGAAGAC", "sample_one_CGTCCTAGGA",
  "sample_one_CGTCGGCAAT", "sample_one_CGTGAGTTAT", "sample_one_CGTGTCAAGC",
  "sample_one_CTAACTTCAG", "sample_one_CTAATAGCGT", "sample_one_CTACACCAGG",
  "sample_one_CTAGCACAAT", "sample_one_CTATGAACGG", "sample_one_CTCAAGGACC",
  "sample_one_CTCACCTGTC", "sample_one_CTCCTATTGT", "sample_one_CTCGCAACGT",
  "sample_one_CTCGTGCCTA", "sample_one_CTGGATTGAC", "sample_one_CTGTAGTCAG",
  "sample_one_CTGTCGCTTC", "sample_one_CTGTCTGTGT", "sample_one_CTTCATATCG",
  "sample_one_CTTGCTGACG", "sample_one_GAAGGATTAG", "sample_one_GAATCGAGCC",
  "sample_one_GACCATCTAA", "sample_one_GACGACCACA", "sample_one_GAGACATCTT",
  "sample_one_GAGCGAGTCA", "sample_one_GAGTAGACCA", "sample_one_GATACGCTTA",
  "sample_one_GATAGACTGT", "sample_one_GATAGAGGCG", "sample_one_GATAGGTCAA",
  "sample_one_GATATCAGGA", "sample_one_GATCTCATTC", "sample_one_GATCTGGTCG",
  "sample_one_GATGAGTGAC", "sample_one_GATGGATACA", "sample_one_GATGTGACAG",
  "sample_one_GATTAAGTCC", "sample_one_GATTGCACGC", "sample_one_GCAAGCGAAT",
  "sample_one_GCAATGTAAG", "sample_one_GCACACTATA", "sample_one_GCACTCGGAA",
  "sample_one_GCACTGCGTT", "sample_one_GCACTTAATC", "sample_one_GCAGGAGATG",
  "sample_one_GCAGTACTGG", "sample_one_GCATATGAGT", "sample_one_GCATCCGATC",
  "sample_one_GCCAAGTACA", "sample_one_GCCACGATTC", "sample_one_GCCATAGGTT",
  "sample_one_GCCATATCGA", "sample_one_GCCGTCAATA", "sample_one_GCCTGGACAT",
  "sample_one_GCGTAATTAC", "sample_one_GCTATTATCC", "sample_one_GCTCAGTAAT",
  "sample_one_GCTGCTTATA", "sample_one_GGAATAAGCA", "sample_one_GGACGATGCT",
  "sample_one_GGCATCGTGA", "sample_one_GGCATTATTG", "sample_one_GGCCGAGATT",
  "sample_one_GGCGCTATAA", "sample_one_GGCGTTAAGT", "sample_one_GGCTATTGAT",
  "sample_one_GGCTGCTACT", "sample_one_GGTAATGTGT", "sample_one_GGTGGTTGGA",
  "sample_one_GGTGTTCACC", "sample_one_GGTTAGATCT", "sample_one_GGTTATGGCG",
  "sample_one_GGTTCACTGG", "sample_one_GGTTGTGCAA", "sample_one_GTAACCAGTA",
  "sample_one_GTAACCTTGG", "sample_one_GTAAGAACCT", "sample_one_GTAAGGCTCC",
  "sample_one_GTAATCCACG", "sample_one_GTATTGTGGA", "sample_one_GTCCGCATCA",
  "sample_one_GTCCTTCGGT", "sample_one_GTCGCTCTCT", "sample_one_GTCGGTGACA",
  "sample_one_GTCTCGAGTG", "sample_one_GTCTCTTAAG", "sample_one_GTCTTCCGAG",
  "sample_one_GTGACTATAC", "sample_one_GTGGTTAATG", "sample_one_GTGTGCCTGT",
  "sample_one_GTGTGTGTCC", "sample_one_GTTCATTGCC", "sample_one_GTTCCGGTGA",
  "sample_one_GTTCGTCGAA", "sample_one_GTTGAATTGG", "sample_one_GTTGATCCGC",
  "sample_one_GTTGTATGCT", "sample_one_TAACCGTAGC", "sample_one_TAACGTCGAT",
  "sample_one_TAAGGTACGG", "sample_one_TACGGACATA", "sample_one_TACTACCGCC",
  "sample_one_TACTGTCAAG", "sample_one_TAGCGAACGC", "sample_one_TAGCGCCAAC",
  "sample_one_TAGGACGCCT", "sample_one_TAGGTTGCAA", "sample_one_TAGTAGTCTC",
  "sample_one_TAGTCCGCTG", "sample_one_TAGTGGAACT", "sample_one_TATCATGCAG",
  "sample_one_TATCGTTACG", "sample_one_TCAAGTGCAG", "sample_one_TCACAGATAC",
  "sample_one_TCACCGCCTA", "sample_one_TCACGCCACT", "sample_one_TCACGTTGGC",
  "sample_one_TCATTGTCCA", "sample_one_TCCACACTAG", "sample_one_TCCACGGTCA",
  "sample_one_TCCACTCGCT", "sample_one_TCCGACTAAC", "sample_one_TCCGTTATCT",
  "sample_one_TCCTAAGAGA", "sample_one_TCCTCTAGTA", "sample_one_TCGAAGCATT",
  "sample_one_TCGAGAGAGC", "sample_one_TCGCACTTGA", "sample_one_TCGCCTACTG",
  "sample_one_TCGCGTAGCA", "sample_one_TCGGCGTTAA", "sample_one_TCTACATCCG",
  "sample_one_TCTCCACATT", "sample_one_TCTCTCCTAT", "sample_one_TCTTGCTCGG",
  "sample_one_TGAACTAACC", "sample_one_TGAAGAAGGT", "sample_one_TGAGCGTTCC",
  "sample_one_TGAGTACGTA", "sample_one_TGGAATGGAG", "sample_one_TGTCATTCGC",
  "sample_one_TGTGCTTCAG", "sample_one_TGTTCAGGAT", "sample_one_TTACACACGT",
  "sample_one_TTACTGTGAC", "sample_one_TTATAGGAGG", "sample_one_TTATCGCGTT",
  "sample_one_TTATGCCGCG", "sample_one_TTCACGGAAG", "sample_one_TTCAGGAGTA",
  "sample_one_TTCCATCGAG", "sample_one_TTCGAGTGAT", "sample_one_TTCTGTACCT",
  "sample_one_TTGGCAATTC", "sample_one_TTGGCTCCAC", "sample_one_TTGGTAACAG",
  "sample_one_TTGGTCAGTA", "sample_one_TTGTCGGCCA", "sample_one_TTGTGTTCGA"
)
stopifnot(identical(sampleNames(sample_1_result), expected_sample_names))

expected_var_labels <- c(
  "WellBC",
  "WellID",
  "NumberOfMTReads",
  "pctMT",
  "NumberOfERCCReads",
  "pctERCC",
  "NumberOfChromReads",
  "pctChrom",
  "NumberOfInputReads",
  "NumberOfMappedReads",
  "PctMappedReads",
  "NumberOfReadsMappedToMultipleLoci",
  "PectOfReadsMappedToMultipleLoci",
  "NumberOfReadsMappedToTooManyLoci",
  "PectOfReadsMappedToTooManyLoci",
  "NumberOfReadsUnmappedTooManyMismatches",
  "PectOfReadsUnmappedTooManyMismatches",
  "NumberOfReadsUnmappedTooShort",
  "PectOfReadsUnmappedTooShort",
  "NumberOfReadsUnmappedOther",
  "PectOfReadsUnmappedOther",
  "ReadsWithValidBarcodes",
  "SequencingSaturation",
  "Q30BasesInCB.UMI",
  "ReadsMappedToTranscriptome.Unique.MultipeGenes",
  "EstimatedNumberOfCells",
  "FractionOfReadsInCells",
  "MeanReadsPerCell",
  "NumberOfUMIs",
  "NumberOfGenes",
  "NumberOfCountedReads",
  "PoolName"
)
stopifnot(identical(varLabels(sample_1_result), expected_var_labels))

read_mm <- function(mapping_dir) {
  market_matrix_file <- file.path(mapping_dir, "Solo.out",
                                  "Gene", "raw", "matrix.mtx")
  result <- readMM(market_matrix_file)
  feature_file <- file.path(mapping_dir, "Solo.out",
                            "Gene", "raw", "features.tsv")
  features <- read.table(feature_file, sep = "\t", header = FALSE,
                         col.names = c("ID", "Name", "Type"))$ID
  rownames(result) <- gsub("\\.\\d+$", "", features)
  barcodes_file <- file.path(mapping_dir,
                             "Solo.out", "Gene", "raw", "barcodes.tsv")
  if (!file.exists(barcodes_file)) {
    stop(paste0("Expected the 'Solo.out/Gene/raw' directory at ",
                mapping_dir, " to contain a 'barcodes.tsv' file."))
  }
  barcodes <- readLines(barcodes_file)
  if (length(barcodes) != 1) {
    stop(paste0("A single STAR Solo folder should only have ",
                "mapped one (1) barcode, but found '",
                length(barcodes), "'for mapping directory ", mapping_dir))
  }
  colnames(result) <- paste0("sample_one_", barcodes)
  return(result)
}
expected_matrices <- lapply(par$star_output, read_mm)
expected_matrix <- as.matrix(do.call(cbind, expected_matrices))
result_counts <- exprs(sample_1_result)
stopifnot(length(setdiff(colnames(expected_matrix),
                         colnames(exprs(sample_1_result)))) == 0)
stopifnot(length(setdiff(rownames(expected_matrix),
                         rownames(exprs(sample_1_result)))) == 0)
expected_matrix_sorted <- expected_matrix[, colnames(exprs(sample_1_result))]
stopifnot(identical(exprs(sample_1_result), expected_matrix_sorted))