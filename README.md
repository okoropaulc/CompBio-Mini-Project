# CompBio-Mini-Project
Computational Biology (COMP 488) Mini Project

This project is purposefully designed to make the use of some bioinformatics software tools easy to run for people with little or no coding skills and those who find it hard running analysis on the command line. 

# The softwares used here are:
PROKKA
TOPHAT2
BOWTIE
CUFFLINKS

# Python Package and version
Biopython
Python3

# How to run the script
python3 mini_project.py --input_path <folder>
The folder is the directory where you want the analysis results and data to be stored.

# Project Description
More specifically and per this project, the analysis contained herein is targetted at quantifying differential gene expression among four different strains of Uropathogenic E.Coli (UPEC). The four strains of UPEC used here are HM27, HM46, HM65, HM69. The genome assemblies are available at:

HM27:  https://www.ncbi.nlm.nih.gov/nuccore/APNU00000000
HM46:  https://www.ncbi.nlm.nih.gov/nuccore/APNY00000000
HM65:  https://www.ncbi.nlm.nih.gov/nuccore/APNX00000000
HM69:  https://www.ncbi.nlm.nih.gov/nuccore/APNV00000000

The transcriptomes for the UPEC strains are available at:

HM27:  https://www.ncbi.nlm.nih.gov/sra/SRX541301
HM46:  https://www.ncbi.nlm.nih.gov/sra/SRX541306
HM65:  https://www.ncbi.nlm.nih.gov/sra/SRX541312
HM69:  https://www.ncbi.nlm.nih.gov/sra/SRX541316

# The python script mini_project.py will implement this analysis by integrating all the above softwares and carrying out the following workflows:

All output of each step below is written out in a logfile, UPEC.log

# Workflow
Firstly, fetch the genome of the four given strains
Secondly, calculate the number of contigs in each assembly
Thirdly, calculate the total length of base pairs in each assembly
Fourth, annotate each assembly using the prokka software
Fifth, since assembled genome in RefSeq for E. coli K-12 (NC_000913) has 4140 CDS and 89 tRNAs annotated. Compare each annotated assembly to the E.coli k-12 annotation.
Sixth, use bowtie2 to build genome index with the UPEC strains.
Seventh, use tophat2 to map each UPEC strain transcriptome to their respective genome.
Eight, use cufflinks to quantify the gene expression.


