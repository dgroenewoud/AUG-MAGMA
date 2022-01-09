# AUG-MAGMA
A controlled approach for incorporating regulatory interactions into gene-set analyses for GWAS data

When using our approach, we ask that you reference our publication as well:

# Background
We provide an Rscript to execute gene scoring and gene-set analysis with MAGMA with the specific aim of comparing between a baseline SNV-to-gene mapping (that 
is, a minimal SNV-to-gene mapping) and an augmented SNV-to-gene mapping (that is, a SNV-to-gene mapping that is built on-top-of a baseline SNV-to-gene mapping). 
For each mapping, MAGMA (a popular gene-set analysis tool for GWAS data: https://doi.org/10.1371/journal.pcbi.1004219) is used to calculate gene scores based on 
SNV-level associations (from GWAS summary statistics) amongst SNVs mapped to each gene. These gene scores are then fed into a competitive gene-set analysis to 
identify collections of genes that are enriched for phenotype association. Gene scores and gene-set scores can be compared between both mappings to evaluate 
the benefits of augmentation. To limit spurious discoveries, we control for non-specific effects with matched, random augmentation (see publication for details).

# Requirements
- Access to a linux server
- An installation of R version 3.5.3 or higher
- An installation of the data.table package (version 1.12.2 or higher), the stringi package (version 1.4.3 or higher), the igraph package (version 1.2.4.1 
  or higher), the foreach package (version 1.4.4 or higher), the parallel package (version 3.5.3 or higher), and the doParallel package (version 1.0.15 or 
  higher)
- An installation of MAGMA (refer to: https://ctg.cncr.nl/software/magma)
- A data set of GWAS summary statistics containing:
  (i) an identifier column with rs identifiers
  (ii) a p-value column
  (iii) a sample-size column (if not provided, set this manually to the study sample-size for every entry)
- A relevant data set of binary files (refer to: https://ctg.cncr.nl/software/magma)
  (i) choose the appropriate population with respect to the GWAS study being analyzed
  (ii) it is possible to use custom binary files with MAGMA at your own risk
- A gene locations file for building SNV-to-gene mappings (refer to: https://ctg.cncr.nl/software/magma):
  (i) ensure that the genome build matches the binary files
  (ii) it is possible to use a custom gene locations file with MAGMA at your own risk
  (iii) it is also possible to avoid this file completely and to build your own SNV-to-gene mappings directly (however, we currently do not implement this)
- A gene-set file (for examples, refer to: http://www.gsea-msigdb.org/gsea/downloads.jsp)

# Prerequisites
Users should be familiar with MAGMA and we recommend that they refer to the MAGMA manual for additional details. The manual can be downloaded from
the MAGMA website (https://ctg.cncr.nl/software/magma).

# Installation
Simply save the Rscript (AUG-MAGMA.R) into any directory of your choice.<br />

Create an empty working directory and, within it, an empty subdirectory called "output". For convenience, we recommend that you also create a subdirectory called
"input" (within the working directory, that is), and that you create five subdirectories within it, namely: "annotations" (for gene locations and SNV-to-gene 
mappings), "binaries" (for a set of binary files), "miscellaneous" (for an optional file containing a list of genes to exclude from the analyses), "sets" (for
a gene-set file), and "sumstats" (for a GWAS summary statistics file). 
  
./workdir/<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;input/<br />
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;annotations/<br />
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;binaries/<br />
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;miscellaneous/<br />
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sets/<br />
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sumstats/<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;output/<br />

# Tutorial
Stage 1. We recommend this tutorial for first-time users to confirm that everything is running as it should. Let's examine the effect of incorporating 10kb 
flanks on-top-of gene bodies in the context of a gene-set analysis for coronary-artery disease GWAS summary statistics. Create the directory structure as described 
above (under the installation heading). Name the working directory "tutorial".
- "Download" and "unzip" the gene locations file from the MAGMA website for human genome build GRCh37 into ./tutorial/input/annotations/: https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip
- "Download" and "unzip" the relevant set of binary files (note, this tutorial concerns a GWAS study performed in a European population) into ./tutorial/input/binaries/: https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
- "Download" a gene-set file (note: since in this tutorial our gene locations file is using entrez identifiers, we should ensure that our gene-set files also uses entrez identifiers) into ./tutorial/input/sets/: http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/c5.go.bp.v7.4.entrez.gmt
- "Download" and "gunzip" GWAS summary statistics for coronary-artery disease into ./tutorial/input/sumstats/: http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz 
- "Rename" the GWAS summary statistics file (e.g. cartery-disease.txt). The file name cannot contain a full-stop. The full-stop is reserved as a delimiter for between the
file name and the file suffix. 

Stage 2. Build two SNV-to-gene mappings. The baseline SNV-to-gene mapping assigns SNVs to genes based on overlap with gene bodies. The augmented SNV-to-gene
mapping assigns SNVs to genes based on overlap with either gene bodies or 10kb flanks. Our publication describes how to build custom SNV-to-gene mappings 
such as mappings that incorporate regulatory interactions. However, for this tutorial we can use the following general command to first build the baseline
SNV-to-gene mapping and then the augmented SNV-to-gene mapping.

</path/to/magma-executable>&nbsp;&nbsp;&nbsp;# define path to magma executable<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--annotate window=X,Y&nbsp;&nbsp;&nbsp;# define upstream (X) and downstream (Y) flank-size in kb<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--snp-loc </path/to/binaries/prefix.bim>&nbsp;&nbsp;&nbsp;# define path to .bim file (one of the binary files)<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--gene-loc </path/to/annotations/gene.loc.file>&nbsp;&nbsp;&nbsp;# define path to gene locations file<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--out </path/to/annotations/genes_uXdY>&nbsp;&nbsp;&nbsp;# define path and name of output file<br />

Stage 3. Run the Rscript from the linux command-line (see below).

# Running AUG-MAGMA

To run AUG-MAGMA, use the following general command in the command line:

nohup<br />
</path/to/R-interpreter>&nbsp;&nbsp;&nbsp;# define path to R interpreter<br />
</path/to/Rscript>&nbsp;&nbsp;&nbsp;# define path to Rscript for execution<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--magma </path/to/magma-executable>&nbsp;&nbsp;&nbsp;# define path to magma executable<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--sumstat </path/to/summary-statistics-file>&nbsp;&nbsp;&nbsp;# define path to summary statistics file<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--sumstat-id Q&nbsp;&nbsp;&nbsp;# Q is column index of column containing rs-identifier<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--sumstat-pval W&nbsp;&nbsp;&nbsp;# W is column index of column containing p-value<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--sumstat-nsample R&nbsp;&nbsp;&nbsp;# R is column index of column containing sample size<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--binaries </path/to/binaries/prefix>&nbsp;&nbsp;&nbsp;# define path to the set of binary files including their common prefix<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--baseline-model </path/to/annotations/baseline-prefix.genes.annot>&nbsp;&nbsp;&nbsp;# define path to baseline SNV-to-gene mapping<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--augmented-model </path/to/annotations/augmented-prefix.genes.annot>&nbsp;&nbsp;&nbsp;# define path to augmented SNV-to-gene mapping<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--gene-set-file </path/to/gene-set-file>&nbsp;&nbsp;&nbsp;# define path to gene-set file<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--output </path/to/output>&nbsp;&nbsp;&nbsp;# define directory for storing output<br />
&<br />

Optional flags *include*:

--cores V&nbsp;&nbsp;&nbsp;# set number of cores V manually (not recommended to exceed default | default is a quarter of total available cores)<br />
--gene-scoring-model top&nbsp;&nbsp;&nbsp;# change the way MAGMA calculates gene scores (not recommended | default is SNP-Wise Mean)<br /> 
--gene-set-format col=A,B&nbsp;&nbsp;&nbsp;# use alt. gene-set format (see MAGMA manual | A/B are index of gene/set column | default is row-based)<br />
--ignore-genes </path/to/gene-list-file>&nbsp;&nbsp;&nbsp;# define path list of genes to exclude from analyses (see MAGMA manual | default is none)<br />
--permutations P&nbsp;&nbsp;&nbsp;# set number of permutations P manually (default is 20)<br />

# Output

Run time varies but will usually not exceed 12 hours. Progress can be monitored in the nohup.out file. All output is stored in the specified output directory
under a subdirectory named (1) after the summary statistics file and then (2) an additional subdirectory named after the augmented SNV-to-gene mapping file. At 
the bottom of this directory structure are three subdirectories: 

The "annotation" subdirectory contains a table in which each each gene is defined by SNVs mapped to it either exclusively via augmentation (a, that is a for 
augmentation) or via both mappings (b, that is b for baseline, since the SNV is already mapped to the gene via the minimal SNV-to-gene mapping). This file serves 
as a reference for any downstream analyses that the user may wish to perform.

The "scores" subdirectory contains gene scores and gene-set scores according to the baseline SNV-to-gene mapping with unpermuted summary statistics ("baseline"), 
the augmented SNV-to-gene mapping with unpermuted summary statistics ("augmented"), and the augmented SNV-to-gene mapping with permuted (EPVP) summary statistics 
("random"). The intermediate output stored within the batches subdirectory (under the random subdirectory) can be ignored. The file names are self explanatory. 
Suffixes:<br />
&nbsp;- Unadjusted gene scores (unadjusted.genes.raw)<br />
&nbsp;- Adjusted gene scores (.adjusted.gsa.genes.out)<br />
&nbsp;- Gene-set scores from competitive gene-set analysis with adjusted gene scores (.adjusted.gsa.out)<br /> 












