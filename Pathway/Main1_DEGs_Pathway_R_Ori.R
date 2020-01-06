
# This code is used to obtain the initial DEGs for control (normal) groups (denoted using All (i.e ctrl) in this procedure) up- and down-upregulated, and normal groups, respectively.

setwd("Your current path")
rm(list = ls())

library(openxlsx)

#################### Step1: data preparation
### file path for counts stage data
Counts_Stage_path <- "~/expression/workspace_one/express_COUNT_Stage_new"

### file path for FPKM stage data
# FPKM_Stage_path <- "~/expression/workspace_one/express_COUNT_Stage_new"
FPKM_Stage_path <- "~/expression/workspace_one/express_FPKM_Stage_new"
### The path of the original cancer staging program
path_ori_program <- "~/Rworkspaces/Cancer_Staging_Project/revision/Cancer_Staging_SSN"
pvalue_cutoff <- 0.05
logFC_cutoff <- 2.5
# logFC_cutoff <- 3 # for LUAD 
# logFC_cutoff <- 2 # for THCA
#################### Step 2: Filter the edgeR results for Main1_DEGs_Stage.R
file_name <- list.files(Counts_Stage_path,pattern = "RData$")
projects  <- sub("(.*)\\..*$", "\\1", file_name)

DEGs_classifier_list_summary <- list()
project_count <- 0
for (project in projects[2:9]) {
  project_count <- project_count + 1
  
  ### load the fpkm stage data for project
  load(paste0(FPKM_Stage_path,'/',project,'.RData'))
  
  #################### Step 2.1: get the filtering DEGs for each stage
  ### load the edgeR results before
  load(paste0(path_ori_program,"/DEGs/",project,"_edgeR_Stages.Rdata"))
  ### 1) edgeR_Stage_I
  index_logFC <- which(abs(edgeR_Stage_I$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_I$p_value < pvalue_cutoff)
  # DEGs_Stage_I <- edgeR_Stage_I$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_I <- edgeR_Stage_I[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_I$p_value,method = "BH")
  DEGs_Stage_I <- cbind(DEGs_Stage_I,adj_pvalue)
  # DEGs_Stage_I_upregulated <- DEGs_Stage_I$symbol[DEGs_Stage_I$log2_fold_change > 0 & DEGs_Stage_I$adj_pvalue < 0.05]
  DEGs_Stage_I_All <- DEGs_Stage_I$symbol[DEGs_Stage_I$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_I_upregulated <- DEGs_Stage_I$symbol[DEGs_Stage_I$log2_fold_change > 0 & DEGs_Stage_I$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_I_downregulated <- DEGs_Stage_I$symbol[DEGs_Stage_I$log2_fold_change < 0 & DEGs_Stage_I$adj_pvalue < pvalue_cutoff]
  cat("The number of up and down regulated DEGs for stage i is:",length(DEGs_Stage_I_upregulated),"and",length(DEGs_Stage_I_downregulated),"\n")
  
  ### 2) edgeR_Stage_II
  index_logFC <- which(abs(edgeR_Stage_II$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_II$p_value < pvalue_cutoff)
  # DEGs_Stage_II <- edgeR_Stage_II$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_II <- edgeR_Stage_II[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_II$p_value,method = "BH")
  DEGs_Stage_II <- cbind(DEGs_Stage_II,adj_pvalue)
  # DEGs_Stage_II_upregulated <- DEGs_Stage_II$symbol[DEGs_Stage_II$log2_fold_change > 0 & DEGs_Stage_II$adj_pvalue < 0.05]
  DEGs_Stage_II_All <- DEGs_Stage_II$symbol[DEGs_Stage_II$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_II_upregulated <- DEGs_Stage_II$symbol[DEGs_Stage_II$log2_fold_change > 0 & DEGs_Stage_II$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_II_downregulated <- DEGs_Stage_II$symbol[DEGs_Stage_II$log2_fold_change < 0 & DEGs_Stage_II$adj_pvalue < pvalue_cutoff]
  cat("The number of up and down regulated DEGs for Stage ii is:",length(DEGs_Stage_II_upregulated),"and",length(DEGs_Stage_II_downregulated),"\n")
  
  ### 1) edgeR_Stage_III
  index_logFC <- which(abs(edgeR_Stage_III$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_III$p_value < pvalue_cutoff)
  # DEGs_Stage_III <- edgeR_Stage_III$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_III <- edgeR_Stage_III[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_III$p_value,method = "BH")
  DEGs_Stage_III <- cbind(DEGs_Stage_III,adj_pvalue)
  # DEGs_Stage_III_upregulated <- DEGs_Stage_III$symbol[DEGs_Stage_III$log2_fold_change > 0 & DEGs_Stage_III$adj_pvalue < 0.05]  
  DEGs_Stage_III_All <- DEGs_Stage_III$symbol[DEGs_Stage_III$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_III_upregulated <- DEGs_Stage_III$symbol[DEGs_Stage_III$log2_fold_change > 0 & DEGs_Stage_III$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_III_downregulated <- DEGs_Stage_III$symbol[DEGs_Stage_III$log2_fold_change < 0 & DEGs_Stage_III$adj_pvalue < pvalue_cutoff]
  cat("The number of up and down regulated DEGs for Stage iii is:",length(DEGs_Stage_III_upregulated),"and",length(DEGs_Stage_III_downregulated),"\n")
  
  ### 1) edgeR_Stage_IV
  index_logFC <- which(abs(edgeR_Stage_IV$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_IV$p_value < pvalue_cutoff)
  # DEGs_Stage_IV <- edgeR_Stage_IV$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_IV <- edgeR_Stage_IV[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_IV$p_value,method = "BH")
  DEGs_Stage_IV <- cbind(DEGs_Stage_IV,adj_pvalue)
  # DEGs_Stage_IV_upregulated <- DEGs_Stage_IV$symbol[DEGs_Stage_IV$log2_fold_change > 0 & DEGs_Stage_IV$adj_pvalue < 0.05]
  DEGs_Stage_IV_All <- DEGs_Stage_IV$symbol[DEGs_Stage_IV$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_IV_upregulated <- DEGs_Stage_IV$symbol[DEGs_Stage_IV$log2_fold_change > 0 & DEGs_Stage_IV$adj_pvalue < pvalue_cutoff]
  DEGs_Stage_IV_downregulated <- DEGs_Stage_IV$symbol[DEGs_Stage_IV$log2_fold_change < 0 & DEGs_Stage_IV$adj_pvalue < pvalue_cutoff]
  cat("The number of up and down regulated DEGs for Stage iv is:",length(DEGs_Stage_IV_upregulated),"and",length(DEGs_Stage_IV_downregulated),"\n")
  
  
  #################### Step 2.2: Integrate the final results
  ### First, get all the genes for all stages
  # DEGs_Name_I <- DEGs_Stage_I
  # DEGs_Name_II <- DEGs_Stage_II
  # DEGs_Name_III <- DEGs_Stage_III
  # DEGs_Name_IV <- DEGs_Stage_IV
  
  ### Then, save the final results
  # save(DEGs_Name_I,DEGs_Name_II,DEGs_Name_III,DEGs_Name_IV,
  #      file = paste0("./DEGs_Ori/",project,"_DEGs_Stages_df.Rdata"))
  save(DEGs_Stage_I_All,DEGs_Stage_II_All,DEGs_Stage_III_All,DEGs_Stage_IV_All,
       DEGs_Stage_I_upregulated,DEGs_Stage_II_upregulated,DEGs_Stage_III_upregulated,DEGs_Stage_IV_upregulated,
       DEGs_Stage_I_downregulated,DEGs_Stage_II_downregulated,DEGs_Stage_III_downregulated,DEGs_Stage_IV_downregulated,
       file = paste0("./DEGs_Ori/",project,"_DEGs_Stages_df.Rdata"))
  cat("The DEGs filtering for pathwway enrichment has been completed!\n")
}

cat('The DEGs save work for pathway enrichment has been completed!')


