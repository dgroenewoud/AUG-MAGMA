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
Stage 1. We recommend the following tutorial for first-time use in order to confirm that everything is running as it should. Let us examine the effect of 
incorporating 10kb flanks on-top-of gene bodies in the context of a gene-set analysis for coronary-artery disease GWAS summary statistics. Create the 
directory structure as described above (second bullet point under the installation heading). Name the working directory "tutorial".
- Download and unzip the gene locations file from the MAGMA website for human genome build GRCh37 into ./tutorial/input/annotations/: https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip
- Download and unzip the relevant set of binary files (note, this tutorial concerns a GWAS study performed in a European population) into ./tutorial/input/binaries/: https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
- Download a gene-set file (note, since in this case the gene locations file uses entrez identifiers, we should use a gene-set file that likewise uses entrez identifiers) into ./tutorial/input/sets/: http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/c5.go.bp.v7.4.entrez.gmt
- Download and gunzip GWAS summary statistics for coronary-artery disease into ./tutorial/input/sumstats/: http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz (note: rename the summary statistics file to cartery_disease.txt)

Stage 2. We first build two SNV-to-gene mappings. The baseline SNV-to-gene mapping assigns SNVs to genes based on overlap with gene bodies only. The 
augmented SNV-to-gene mapping assigns SNVs to genes based on overlap with gene bodies or 10kb flanks. Refer to our publication for instructions on how 
to build custom SNV-to-gene mappings such as mappings that incorporate regulatory interactions. 

For this tutorial, we can simply use the following command to first build the baseline SNV-to-gene mapping and then to build the augmented SNV-to-gene 
mapping (note: X and Y represent the upstream and downstream flank size in kb, respectively):

</path/to/magma-executable>&nbsp;&nbsp;&nbsp;# define path to magma executable<br />
--annotate window=X,Y&nbsp;&nbsp;&nbsp;# define upstream (X) and downstream (Y) flank-size in kb<br />
--snp-loc </path/to/binaries/prefix.bim>&nbsp;&nbsp;&nbsp;# define path to .bim file (one of the binary files)<br />
--gene-loc </path/to/annotations/gene.loc.file>&nbsp;&nbsp;&nbsp;# define path to gene locations file<br />
--out </path/to/annotations/genes_uXdY>&nbsp;&nbsp;&nbsp;# define path and name of output file<br />

Stage 3. We are now ready to run the Rscript. The entire process usually takes no more than 12 hours to complete. 

We execute the Rscript from the linux command line. The basic command is:

nohup<br />
</path/to/R-interpreter>&nbsp;&nbsp;&nbsp;# define path to R interpreter<br />
</path/to/Rscript>&nbsp;&nbsp;&nbsp;# define path to Rscript for execution<br />
--magma </path/to/magma-executable>&nbsp;&nbsp;&nbsp;# define path to magma executable<br />
--sumstat </path/to/summary-statistics-file>&nbsp;&nbsp;&nbsp;# define path to summary statistics file<br />
--sumstat-id Q&nbsp;&nbsp;&nbsp;# Q is column index of column containing rs-identifier<br />
--sumstat-pval W&nbsp;&nbsp;&nbsp;# W is column index of column containing p-value<br />
--sumstat-nsample R&nbsp;&nbsp;&nbsp;# R is column index of column containing sample size<br />
--binaries </path/to/binaries/prefix>&nbsp;&nbsp;&nbsp;# define path to binary files including common prefix<br />
--baseline-model </path/to/annotations/baseline-prefix.genes.annot>&nbsp;&nbsp;&nbsp;# define path to baseline SNV-to-gene mapping<br />
--augmented-model </path/to/annotations/augmented-prefix.genes.annot>&nbsp;&nbsp;&nbsp;# define path to augmented SNV-to-gene mapping<br />
--gene-set-file </path/to/gene-set-file>&nbsp;&nbsp;&nbsp;# define path to gene-set file<br />
--output </path/to/output>&nbsp;&nbsp;&nbsp;# define directory for storing output<br />
&<br />

Optional flags include:

--cores V&nbsp;&nbsp;&nbsp;# set number of cores V manually (not recommended)<br />
--gene-scoring-model top&nbsp;&nbsp;&nbsp;# this changes the way MAGMA calculates gene scores (not recommended)<br /> 
--gene-set-format col=A,B&nbsp;&nbsp;&nbsp;# A is index of gene column and B is index of set column (alternative gene-set file format)<br />
--ignore-genes </path/to/gene-list-file>&nbsp;&nbsp;&nbsp;# define path to file containing list of genes to exclude from analyses<br />
--permutations P&nbsp;&nbsp;&nbsp;# set number of permutations P manually<br />











