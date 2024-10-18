library(waldo)
library(Biobase)

old <- readRDS('/home/di/code/htrnaseq/eset_Exp_AGGATAGC.rds')
new <- readRDS('/home/di/code/htrnaseq/output_old_comparison/eset.Exp_AGGATAGC.rds')

# Sort sample names
new <- new[, sampleNames(old)]
var_labels_to_ignore <- c("Compound", "ProjectName",
                          "Remark", "Dose", "DoseUnit", "ELN_Number",
                          "IDPdescription", "IDT", "ExperimentID",
                          "Location", "index", "Quadrant", "RunID",
                          "SampleName", "ReplicateID", "Species",
                          "StudyName", "StudyType", "Subject_ID",
                          "UniqueSampleID", "WellID", "WellID_96well",
                          "WellType", "row_num", "position",
                          "LQ_PC1", "LQ_PC2", "Imp_LQ_PC1", "Imp_LQ_PC2",
                          "general_outlier", "LQ_outlier", "row", "control",
                          "LC_outlier", "column", "qcPassFail",
                          "qcScore_ERCC_hardcutoff",
                          "qcScore_MappedReadsPerUMIs_datadriven",
                          "qcScore_MappedReadsPerUMIs_hardcutoff",
                          "qcScore_MT_datadriven",
                          "qcScore_MT_hardcutoff",
                          "qcScore_NGenes_datadriven",
                          "qcScore_NGenes_hardcutoff",
                          "qcScore_NUMIs_datadriven",
                          "qcScore_NUMIs_hardcutoff",
                          "qcScore_RunBasedNGenesDipTest",
                          "qcTotalScore",
                          # The labels used to be part of the quality_control,
                          # but they could potentionally still be calculated
                          "nMappedReads_per_UMI",
                          "nUMI_per_Gene",
                          "pctCountedReads")
var_labels_to_keep <- varLabels(old)[!varLabels(old) %in% var_labels_to_ignore]
var_labels_to_compare <- c(varLabels(new), var_labels_to_ignore)
phenoData(old) <- phenoData(old)[, varLabels(new)]
f_data_to_ignore = c("ENTREZID", "GENENAME")
f_data_to_select = colnames(fData(old))[!colnames(fData(old)) %in% f_data_to_ignore]
fData(old) <- fData(old)[f_data_to_select]
fData(new) <- fData(new)[colnames(fData(old))]
print(compare(old, new, tolerance=1e-6))
