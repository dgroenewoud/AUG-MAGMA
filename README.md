# EPM-MAGMA
Incorporating regulatory interactions into gene-set analyses for GWAS data

We provide an Rscript to execute gene scoring and gene-set analysis with MAGMA with the aim of comparing between a baseline SNV-to-gene mapping and
an augmented SNV-to-gene mapping. 

A baseline SNV-to-gene mapping is a minimal SNV-to-gene mapping. Augmentation refers to the assignment of additional SNVs to genes. For example, a 
baseline SNV-to-gene mapping may be built by assigning SNVs to genes based on proximity (e.g., assign SNVs to a gene if they overlap with the gene 
body). An augmented SNV-to-gene mapping in this case may be built from adding flanks to genes (e.g., assign any SNVs to a gene if they are located 
within either 10kb upstream of the transcription start-site or 10kb downstream of the transcription end-site). MAGMA, a popular gene-set analysis
tool for GWAS data, is then used to execute gene scoring and gene-set analysis separately with each SNV-to-gene mapping. Gene scores and gene-set 
scores are computed based on the SNV-level associations amongst SNVs mapped to each gene. Scores can ultimately be compared between both mappings.

We demonstrated (see publication) that it is important to control such comparative analyses with matched, random augmentation. Here, we implement 
such a control strategy, which we call extragenic p-value push (EPVP). By default, 20 permutations of EPVP are executed wherein the extragenic SNVs 
mapped to each gene (that is, SNVs mapped to a gene exclusively via augmentation) are assigned background SNV-level associations (that is, SNV-level 
associations sampled from random locations in the genome). This is analogous to sampling matched, random augmentation from the genome whilst striving 
to maintain potential confounders of gene scores and gene set scores (such as the number of SNVs mapped to each gene). As before, gene scores and 
gene-set scores are then calculated but under these permuted conditions. One can then evaluate whether or not a improvement in the association of a 
gene score or gene-set score with the genuine augmentation (that is, relative to the baseline) is specific to it (that is, whether or not a comparable 
improvement occurs with matched, random augmentation). For details, refer to our publication.   






