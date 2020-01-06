
# This code is used to filter the edgeR results

setwd("Your current path")
rm(list = ls())

#################### Step1: data preparation
### file path for counts stage data
Counts_Stage_path <- "~/expression/workspace_one/express_COUNT_Stage_new"

pvalue_cutoff <- 0.05
logFC_cutoff <- 6
#################### Step 2: Filter the edgeR results for Main1_DEGs_Stage.R
file_name <- list.files(Counts_Stage_path,pattern = "RData$")
projects  <- sub("(.*)\\..*$", "\\1", file_name)

for (project in projects[2:9]) {
  ### load the edgeR results before
  load(paste0("./DEGs/",project,"_edgeR_Stages.Rdata"))
  ### 1) edgeR_Stage_I
  index_logFC <- which(abs(edgeR_Stage_I$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_I$p_value < pvalue_cutoff)
  # DEGs_Stage_I <- edgeR_Stage_I$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_I <- edgeR_Stage_I[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_I$p_value,method = "BH")
  DEGs_Stage_I <- cbind(DEGs_Stage_I,adj_pvalue)
  
  # DEGs_Stage_I_upregulated <- DEGs_Stage_I$symbol[DEGs_Stage_I$log2_fold_change > 0 & DEGs_Stage_I$adj_pvalue < 0.05]
  DEGs_Stage_I <- DEGs_Stage_I$symbol[DEGs_Stage_I$adj_pvalue < pvalue_cutoff]
  
  ### 2) edgeR_Stage_II
  index_logFC <- which(abs(edgeR_Stage_II$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_II$p_value < pvalue_cutoff)
  # DEGs_Stage_II <- edgeR_Stage_II$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_II <- edgeR_Stage_II[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_II$p_value,method = "BH")
  DEGs_Stage_II <- cbind(DEGs_Stage_II,adj_pvalue)
  
  # DEGs_Stage_II_upregulated <- DEGs_Stage_II$symbol[DEGs_Stage_II$log2_fold_change > 0 & DEGs_Stage_II$adj_pvalue < 0.05]
  DEGs_Stage_II <- DEGs_Stage_II$symbol[DEGs_Stage_II$adj_pvalue < pvalue_cutoff]
  
  ### 1) edgeR_Stage_III
  index_logFC <- which(abs(edgeR_Stage_III$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_III$p_value < pvalue_cutoff)
  # DEGs_Stage_III <- edgeR_Stage_III$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_III <- edgeR_Stage_III[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_III$p_value,method = "BH")
  DEGs_Stage_III <- cbind(DEGs_Stage_III,adj_pvalue)
  
  # DEGs_Stage_III_upregulated <- DEGs_Stage_III$symbol[DEGs_Stage_III$log2_fold_change > 0 & DEGs_Stage_III$adj_pvalue < 0.05]
  DEGs_Stage_III <- DEGs_Stage_III$symbol[DEGs_Stage_III$adj_pvalue < pvalue_cutoff]
  
  ### 1) edgeR_Stage_IV
  index_logFC <- which(abs(edgeR_Stage_IV$log2_fold_change)> logFC_cutoff)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=pvalue_cutoff)
  index_Pvalue <- which(edgeR_Stage_IV$p_value < pvalue_cutoff)
  # DEGs_Stage_IV <- edgeR_Stage_IV$symbol[intersect(index_logFC,index_Pvalue)]
  DEGs_Stage_IV <- edgeR_Stage_IV[intersect(index_logFC,index_Pvalue),]
  adj_pvalue <- p.adjust(DEGs_Stage_IV$p_value,method = "BH")
  DEGs_Stage_IV <- cbind(DEGs_Stage_IV,adj_pvalue)
  
  # DEGs_Stage_IV_upregulated <- DEGs_Stage_IV$symbol[DEGs_Stage_IV$log2_fold_change > 0 & DEGs_Stage_IV$adj_pvalue < 0.05]
  DEGs_Stage_IV <- DEGs_Stage_IV$symbol[DEGs_Stage_IV$adj_pvalue < pvalue_cutoff]
  
  save(DEGs_Stage_I,DEGs_Stage_II,DEGs_Stage_III,DEGs_Stage_IV,
       file = paste0("./DEGs/",project,"_DEGs_Stages_df.Rdata"))
  
  cat("The edgeR analysis has been completed!\n")
  
}