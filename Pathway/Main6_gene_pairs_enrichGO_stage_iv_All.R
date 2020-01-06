# This code is used to conduct the gene pairs enrichments
setwd("Your current path")
rm(list = ls())
library(openxlsx)

#################### Step 1: load the data
projects <- c("BLCA_Stage","BRCA_Stage","COAD_Stage","HNSC_Stage","KIRC_Stage","KIRP_Stage","LUAD_Stage","STAD_Stage","THCA_Stage")

PCC_cutoff=c(0.60, 0.70, 0.80)
for(pcc_cutoff_item in PCC_cutoff){
  cat("The work for pcc_cutoff_item ",pcc_cutoff_item,"Starts!\n")
  Genepairs_xlsx_ctrl_result <- list()
  Genepairs_xlsx_experi_up_result <- list()
  Genepairs_xlsx_experi_down_result <- list()
  
  projects_count <- 0
  for(project in projects[c(2:9)]){
    projects_count <- projects_count + 1
    ### load the enrihed pathways
    load(paste("./main5_enrichGO_results_Initial_All/",project,"_enrichGO_results.Rdata",sep = ""))
    
    ### load the differential structure for control group and experimental group across each stage
    load(paste0("./CEGs_Part2/",pcc_cutoff_item,"/",project,"_CEGs_Stages_Matrix.Rdata"))
    
    ###################################### Part A: For Stage I
    #################### Step 2:construct reference information for control group of stage i
    # the initial structure
    DEGs_len_ctrl <- dim(CEG_Matrix_N_IV)[1]
    
    ### the matrix form
    # CEG_Matrix_N_IV
    # CEG_Matrix_Stage_IV
    ### the vector form
    ctrl_IV <- which(CEG_Matrix_N_IV == 100)
    
    ### the coming result step: the length of the non zero term is :
    N_ctrl <- DEGs_len_ctrl * (DEGs_len_ctrl)/2
    K_ctrl <- length(ctrl_IV)
    # NK_ctrl <- c(N_ctrl,K_ctrl) ################################################## backgroud
    NK_ctrl <- paste(K_ctrl,"/",N_ctrl,sep = "") ################################################## backgroud
    
    #################### Step 3: construct the experimetal information using the enrichGO_result_IV.result daat
    GO_pvalue_ctrl <- list()
    nk_ctrl <- list()
    GOcount_ctrl <- list()
    padjust_ctrl <- NULL
    count_ctrl <- 0
    for(qq in 1: dim(enrichGO_result_IV.result)[1]){
      count_ctrl <- count_ctrl + 1
      ### get the corresponding geneID for the term qq
      geneID_each <- unlist(strsplit(enrichGO_result_IV.result$geneID[qq] ,split="/"))
      ### get sub reference index matrix for the ctrl data: CEG_Matrix_N_IV
      CEG_Matrix_N_IV_each <- CEG_Matrix_N_IV[geneID_each,geneID_each]
      len_each <- length(geneID_each)
      
      n_ctrl <- len_each*(len_each+1)/2
      # k_ctrl <- length(which(is.na(CEG_Matrix_N_IV_each)==FALSE))
      k_ctrl <-length(which(CEG_Matrix_N_IV_each >=100))
      # M <- K
      # m+n <- N
      # x <- k
      # k <- n
      # ctrl_pvalue <- dhyper(k_ctrl, K_ctrl, N_ctrl-K_ctrl, n_ctrl)
      ctrl_pvalue <- 1 - phyper(k_ctrl-1, K_ctrl, N_ctrl-K_ctrl, n_ctrl)
      # nk_ctrl[[count_ctrl]] <- c(n_ctrl,k_ctrl)################################################## gene pairs ratio
      nk_ctrl[[count_ctrl]] <- paste(k_ctrl,"/",n_ctrl,sep = "")################################################## gene pairs ratio
      GOcount_ctrl[[count_ctrl]] <-  k_ctrl
      GO_pvalue_ctrl[[count_ctrl]] <- ctrl_pvalue ################################################## pvalue
      
    }
    
    GO_pvalue_ctrl <- unlist(GO_pvalue_ctrl) 
    padjust_ctrl <- p.adjust(GO_pvalue_ctrl, method = "BH", n = length(GO_pvalue_ctrl)) ############ p.adjust
    
    ### so we get the following the value
    # nk_ctrl
    # NK_ctrl
    # GO_pvalue_ctrl  : pvalue
    # padjust_ctrl
    # 
    ######  we should construct a dataframe for the new enrichment results of gene pairs
    goterm_len <- dim(enrichGO_result_IV.result)[1]
    GPRatio_ctrl <- unlist(nk_ctrl)
    BgRatio_ctrl <- rep(NK_ctrl,goterm_len)
    # padjust_ctrl
    GOcount_ctrl <- unlist(GOcount_ctrl)
    GenePairs_enrichGO_result_IV.result_ctrl <- NULL
    GenePairs_enrichGO_result_IV.result_ctrl <- data.frame("ID" = enrichGO_result_IV.result$ID,
                                                                  "Description" = enrichGO_result_IV.result$Description,
                                                                  "GPRatio" = GPRatio_ctrl,
                                                                  "BgRatio" = BgRatio_ctrl,
                                                                  "pvalue" = GO_pvalue_ctrl,
                                                                  "p.adjust" = padjust_ctrl,
                                                                  "GOcount" = GOcount_ctrl)
    rownames(GenePairs_enrichGO_result_IV.result_ctrl) <- GenePairs_enrichGO_result_IV.result_ctrl$ID
    
    GO_term_ctrl <- enrichGO_result_IV.result[which(GO_pvalue_ctrl < 0.05),]
    GenePairs_GO_term_ctrl <- GenePairs_enrichGO_result_IV.result_ctrl[which(GO_pvalue_ctrl < 0.05),]
    GenePairs_GO_term_ctrl <- GenePairs_GO_term_ctrl[order(GenePairs_GO_term_ctrl$pvalue),]
    save(GO_term_ctrl,GenePairs_GO_term_ctrl,file = paste("./main6_enrichGO_results_Final/085/",pcc_cutoff_item,"/","stage_iv/",project,"_gene_pairs_GOterms_ctrl_085_stageiv.Rdata",sep = ""))
    
    
    Genepairs_xlsx_ctrl_result[[projects_count]] <- GenePairs_GO_term_ctrl
    cat("This work for ",project," is completed!\n")
  }
  
  names(Genepairs_xlsx_ctrl_result) <- projects[c(2:9)]
  write.xlsx(Genepairs_xlsx_ctrl_result,file = paste("./main6_enrichGO_results_Final/085/",pcc_cutoff_item,"/","stage_iv/","List_gene_pairs_GOterms_ctrl_085_stageiv.xlsx",sep = ""))
  cat("This work for stage i is completed!\n")
}



  



