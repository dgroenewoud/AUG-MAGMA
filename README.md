# EPM-MAGMA
Incorporating regulatory interactions into gene-set analyses for GWAS data

# Background
We provide an Rscript to execute gene scoring and gene-set analysis with MAGMA with the aim of comparing between a baseline SNV-to-gene mapping and
an augmented SNV-to-gene mapping. A baseline SNV-to-gene mapping is a minimal SNV-to-gene mapping (e.g., assign SNVs to a gene if they overlap with
the gene body). An augmented SNV-to-gene mapping is built on-top-of a baseline SNV-to-gene mapping (e.g., assign additional SNVs to a gene if they 
are located either within 10kb upstream of the transcription start-site or 10kb downstream of the transcription end-site). MAGMA (a popular gene-set
analysis tool for GWAS data: https://doi.org/10.1371/journal.pcbi.1004219) is then used to calculate gene scores based on SNV-level associations (from 
GWAS summary statistics) amongst SNVs mapped to each gene. These gene scores are then fed into a competitive gene-set analysis to identify collections
of genes that are enriched for phenotype association. Gene scores and gene-set scores can be compared between both mappings to evaluate the benefits
of augmentation. We demonstrated that, in order to limit spurious discoveries, it is important to control such comparative analyses with matched, 
random augmentation, and we implement such a control (EPVP, see our publication for details) in this Rscript.

# Requirements
- Access to a linux server
- An installation of R version 3.5.3 or higher
- An installation of the data.table package version 1.12.2 or higher
- An installation of MAGMA (refer to: https://ctg.cncr.nl/software/magma)
- A data set of GWAS summary statistics containing:
  (i) an identifier column with rs identifiers
  (ii) a p-value column
  (iii) a sample-size column (if not provided, set this to the study sample-size for every entry)
- A relevant data set of binary files (refer to: https://ctg.cncr.nl/software/magma)
  (i) choose the appropriate population with respect to the GWAS study being analyzed
  (ii) it is possible to use custom binary files
- A gene locations file for building SNV-to-gene mappings (refer to: https://ctg.cncr.nl/software/magma):
  (i) ensure that the genome build matches the binary files
  (ii) it is possible to use a custom gene locations file
  (iii) it is also possible to avoid this file completely and to build your own SNV-to-gene mappings (however, this is not implemented in our Rscript)
- A gene-set file (for examples, refer to: http://www.gsea-msigdb.org/gsea/downloads.jsp)

# Installation
- Simply save the Rscript (EPM-MAGMA.R) into a directory of your choice.
- Create an empty working directory with an "input" and an "output" subdirectory. The input subdirectory should contain five more subdirectories, as follows: 
  "annotations" (for gene locations and SNV-to-gene mappings), "binaries" (for a set of binary files), "miscellaneous" (for an optional text file containing 
  a list of genes to exclude from the analyses), "sets" (for a gene-set file), and "sumstats" (for a GWAS summary statistics file).  

# Tutorial
We recommend the following tutorial for first-time use in order to confirm that everything is running as it should. Let us examine the effect of 
incorporating 10kb flanks on-top-of gene bodies in the context of a gene-set analysis for coronary-artery disease GWAS summary statistics. Create the 
directory structure as described above (under the installation heading). Name the working directory "tutorial".
- Download GWAS summary statistics for coronary-artery disease into ./tutorial/input/sumstats/: http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz
- Download the gene locations file from the MAGMA website for human genome build GRCh37 into ./tutorial/input/annotations/: https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip




