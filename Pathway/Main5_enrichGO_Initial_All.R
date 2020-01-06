
# This code is used to obtain the initial enrichment analysis for control(ctrl) groups.
setwd("Your current path")
rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)

#################### Step 1: load the data
projects <- c("BLCA_Stage","BRCA_Stage","COAD_Stage","HNSC_Stage","KIRC_Stage","KIRP_Stage","LUAD_Stage","STAD_Stage","THCA_Stage")

#################### Step 2: conduct the enrichGO analysis
for(project in projects[c(2:9)]){
  #################### Step 2.1: load the DEGs for project
  load(paste0("./DEGs_Ori/",project,"_DEGs_Stages_df.Rdata"))
  ### conduct enrichGO() function for stage i
  enrichGO_result_I <- enrichGO(gene = DEGs_Stage_I_All, 
                              ont = "BP", 
                              keyType  = 'SYMBOL',
                              OrgDb ="org.Hs.eg.db",
                              pvalueCutoff = 0.05)
  enrichGO_result_I.result <- enrichGO_result_I@result
  
  enrichGO_result_I_sim_085  <- simplify(enrichGO_result_I, 
                                      cutoff=0.85, 
                                      by="p.adjust", 
                                      select_fun=min)
  enrichGO_result_I_sim_085.result <- enrichGO_result_I_sim_085@result
  
  ### conduct enrichGO() function for stage ii
  enrichGO_result_II <- enrichGO(gene = DEGs_Stage_II_All, 
                                ont = "BP", 
                                keyType  = 'SYMBOL',
                                OrgDb ="org.Hs.eg.db",
                                pvalueCutoff = 0.05)
  enrichGO_result_II.result <- enrichGO_result_II@result
  
  enrichGO_result_II_sim_085  <- simplify(enrichGO_result_II, 
                                         cutoff=0.85, 
                                         by="p.adjust", 
                                         select_fun=min)
  enrichGO_result_II_sim_085.result <- enrichGO_result_II_sim_085@result
  
  
  ### conduct enrichGO() function for stage iii
  enrichGO_result_III <- enrichGO(gene = DEGs_Stage_III_All, 
                                 ont = "BP", 
                                 keyType  = 'SYMBOL',
                                 OrgDb ="org.Hs.eg.db",
                                 pvalueCutoff = 0.05)
  enrichGO_result_III.result <- enrichGO_result_III@result
  
  enrichGO_result_III_sim_085  <- simplify(enrichGO_result_III, 
                                          cutoff=0.85, 
                                          by="p.adjust", 
                                          select_fun=min)
  enrichGO_result_III_sim_085.result <- enrichGO_result_III_sim_085@result
  
  
  ### conduct enrichGO() function for stage iv
  enrichGO_result_IV <- enrichGO(gene = DEGs_Stage_IV_All, 
                                  ont = "BP", 
                                  keyType  = 'SYMBOL',
                                  OrgDb ="org.Hs.eg.db",
                                  pvalueCutoff = 0.05)
  enrichGO_result_IV.result <- enrichGO_result_IV@result
  
  enrichGO_result_IV_sim_085  <- simplify(enrichGO_result_IV, 
                                           cutoff=0.85, 
                                           by="p.adjust", 
                                           select_fun=min)
  enrichGO_result_IV_sim_085.result <- enrichGO_result_IV_sim_085@result
  
  
  
  save(enrichGO_result_I.result,enrichGO_result_I_sim_085.result,
       enrichGO_result_II.result,enrichGO_result_II_sim_085.result,
       enrichGO_result_III.result,enrichGO_result_III_sim_085.result,
       enrichGO_result_IV.result,enrichGO_result_IV_sim_085.result,
       file = paste("./main5_enrichGO_results_Initial_All/",project,"_enrichGO_results.Rdata",sep = ""))
  
  cat("This work for ",project," is completed!\n")  
 
  }
  
  cat("This work is completed!\n")  
  
  
