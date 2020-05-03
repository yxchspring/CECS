# CECS
Co-Expression based Cancer Staging and Application

#################### Code description and introduction: ####################

This codes are divided into two parts. 
1) The first part is used for conducting cancer staging classification based on perturbed co-expression networks and,
2) the second part is used for conducting pathway enrichment analysis based on co-expressed genes.


#################### Part1 Classifier:cancer staging classification ####################
1. Main1_DEGs_Stages.R

This code is used to conduct differentially expression analysis (normal VS each stage)

2. Main2_DEGs_Filter.R

This code is used to filter the edgeR results by set the specific cutoff values for pvalue and logFC.

3. Main3_PCC_Stages_CECS.R

This code is used to conduct the co-expression networks(i.e. PCC) and perturbed PCC for each Stage and complete the classification task.


#################### Part2 Pathway: pathway enrichment analysis ####################
1. Main1_DEGs_Pathway_R_Ori.R

This code is used to obtain a) initial DEGs, the b) up- and c) down-regulated DEGs by comparing with each stage and its normal tissues, respectively.

2. Main2_CEGs_g_Up_or_Down_Part1.R
This code is used to obtian the 
a) PCC network using initial DEGs (i.e. Ref_PCC_Stage_N network), 
b) the perturbed PCC using up-regulated DEGs (i.e.  PCC_Stage_up), 
c) the perturbed PCC using down-regulated DEGs (i.e. PCC_Stage_down), 
for each stage, respectively.

3. Main2_CEGs_g_Up_or_Down_Part2.R

This code is used to otain the final significantly differenatilal locations by setting the PCC_cutoff value for each stage, respectively.
Finally, the following results are obtained for each stage, repectively:
CEG_Matrix_N_I, CEG_Matrix_Stage_up, CEG_Matrix_Stage_down.

4. Main5_enrichGO_Initial_All.R
This code is used to obtain the initial enrichment analysis for normal group using a) initial DEGs.

5. Main5_enrichGO_Initial_down.R
This code is used to obtain the initial enrichment analysis for down-regulated group using c) down-regulated DEGs.

6. Main5_enrichGO_Initial_up.R
This code is used to obtain the initial enrichment analysis for up-regulated group using b) up-regulated DEGs.

7.Main6_gene_pairs_enrichGO_stage_i_All.R

This code is used to conduct the pathway enrichment analysis based on co-expressed genes for normal group.

8. Main6_gene_pairs_enrichGO_stage_i_down.R

This code is used to conduct the pathway enrichment analysis based on co-expressed genes for down-regulated group.


9. Main6_gene_pairs_enrichGO_stage_i_up.R
This code is used to conduct the pathway enrichment analysis based on co-expressed genes for up-regulated group.

Note: The operations of the remaining three stages are consistent with those of the stage i.


















