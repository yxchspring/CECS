
# This code is used to conduct the pathway enrichment analysis based on co-expressed genes for down-regulated group.
setwd("Your current path")
rm(list = ls())

library(openxlsx)

#################### Step 1: load the data
projects <- c("BLCA_Stage","BRCA_Stage","COAD_Stage","HNSC_Stage","KIRC_Stage","KIRP_Stage","LUAD_Stage","STAD_Stage","THCA_Stage")

PCC_cutoff=c(0.60, 0.70, 0.80)
for(pcc_cutoff_item in PCC_cutoff){
  cat("The work for pcc_cutoff_item ",pcc_cutoff_item,"Starts!\n")
  Genepairs_xlsx_exper_result <- list()
  
  projects_count <- 0
  for(project in projects[c(2:9)]){
    projects_count <- projects_count + 1
    ### load the enrihed pathways
    load(paste("./main5_enrichGO_results_Initial_down/",project,"_enrichGO_results.Rdata",sep = ""))
    
    ### load the differential structure for control group and experimental group across each stage
    load(paste0("./CEGs_Part2/",pcc_cutoff_item,"/",project,"_CEGs_Stages_Matrix.Rdata"))
    
    ###################################### Part A: For Stage I
    #################### Step 2:construct reference information for control group of stage i
    # the initial structure
    DEGs_len_exper <- dim(CEG_Matrix_Stage_I_down)[1]
    
    ### the matrix form
    # CEG_Matrix_N_I
    # CEG_Matrix_Stage_I
    ### the vector form
    exper_I_down <- which(CEG_Matrix_Stage_I_down == 100)
    
    ### the coming result step: the length of the non zero term is :
    N_exper <- DEGs_len_exper * (DEGs_len_exper)/2
    K_exper <- length(exper_I_down)
    # NK_exper <- c(N_exper,K_exper) ################################################## backgroud
    NK_exper <- paste(K_exper,"/",N_exper,sep = "") ################################################## backgroud
    
    #################### Step 3: construct the experimetal information using the enrichGO_result_I.result daat
    GO_pvalue_exper <- list()
    nk_exper <- list()
    GOcount_exper <- list()
    padjust_exper <- NULL
    count_exper <- 0
    for(qq in 1: dim(enrichGO_result_I.result)[1]){
      count_exper <- count_exper + 1
      ### get the corresponding geneID for the term qq
      geneID_each <- unlist(strsplit(enrichGO_result_I.result$geneID[qq] ,split="/"))
      ### get sub reference index matrix for the exper data: CEG_Matrix_Stage_I_down
      CEG_Matrix_Stage_I_down_each <- CEG_Matrix_Stage_I_down[geneID_each,geneID_each]
      len_each <- length(geneID_each)
      
      n_exper <- len_each*(len_each+1)/2
      # k_exper <- length(which(is.na(CEG_Matrix_Stage_I_down_each)==FALSE))
      k_exper <-length(which(CEG_Matrix_Stage_I_down_each >=100))
      # M <- K
      # m+n <- N
      # x <- k
      # k <- n
      # exper_pvalue <- dhyper(k_exper, K_exper, N_exper-K_exper, n_exper)
      exper_pvalue <- 1 - phyper(k_exper-1, K_exper, N_exper-K_exper, n_exper)
      # nk_exper[[count_exper]] <- c(n_exper,k_exper)################################################## gene pairs ratio
      nk_exper[[count_exper]] <- paste(k_exper,"/",n_exper,sep = "")################################################## gene pairs ratio
      GOcount_exper[[count_exper]] <-  k_exper
      GO_pvalue_exper[[count_exper]] <- exper_pvalue ################################################## pvalue
      
    }
    
    GO_pvalue_exper <- unlist(GO_pvalue_exper) 
    padjust_exper <- p.adjust(GO_pvalue_exper, method = "BH", n = length(GO_pvalue_exper)) ############ p.adjust
    
    ### so we get the following the value
    # nk_exper
    # NK_exper
    # GO_pvalue_exper  : pvalue
    # padjust_exper
    # 
    ######  we should construct a dataframe for the new enrichment results of gene pairs
    goterm_len <- dim(enrichGO_result_I.result)[1]
    GPRatio_exper <- unlist(nk_exper)
    BgRatio_exper <- rep(NK_exper,goterm_len)
    # padjust_exper
    GOcount_exper <- unlist(GOcount_exper)
    GenePairs_enrichGO_result_I.result_exper <- NULL
    GenePairs_enrichGO_result_I.result_exper <- data.frame("ID" = enrichGO_result_I.result$ID,
                                                          "Description" = enrichGO_result_I.result$Description,
                                                          "GPRatio" = GPRatio_exper,
                                                          "BgRatio" = BgRatio_exper,
                                                          "pvalue" = GO_pvalue_exper,
                                                          "p.adjust" = padjust_exper,
                                                          "GOcount" = GOcount_exper)
    rownames(GenePairs_enrichGO_result_I.result_exper) <- GenePairs_enrichGO_result_I.result_exper$ID
    
    GO_term_exper <- enrichGO_result_I.result[which(GO_pvalue_exper < 0.05),]
    GenePairs_GO_term_exper <- GenePairs_enrichGO_result_I.result_exper[which(GO_pvalue_exper < 0.05),]
    GenePairs_GO_term_exper <- GenePairs_GO_term_exper[order(GenePairs_GO_term_exper$pvalue),]
    save(GO_term_exper,GenePairs_GO_term_exper,file = paste("./main6_enrichGO_results_Final/085/",pcc_cutoff_item,"/","stage_i/",project,"_gene_pairs_GOterms_exper_085_stagei_down.Rdata",sep = ""))
    
    
    Genepairs_xlsx_exper_result[[projects_count]] <- GenePairs_GO_term_exper
    cat("This work for ",project," is completed!\n")
  }
  
  names(Genepairs_xlsx_exper_result) <- projects[c(2:9)]
  write.xlsx(Genepairs_xlsx_exper_result,file = paste("./main6_enrichGO_results_Final/085/",pcc_cutoff_item,"/","stage_i/","List_gene_pairs_GOterms_exper_085_stagei_down.xlsx",sep = ""))
  cat("This work for stage i is completed!\n")
}







