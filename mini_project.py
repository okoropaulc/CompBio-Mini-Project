#!/usr/bin/python3

import os
from Bio import SeqIO
from Bio.Seq import Seq

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_path", type = str, action = "store", required = True, help = "Path to folder containing your data")
args = parser.parse_args()

path = args.input_path

#os.chdir("/home/pokoro/Paul_Okoro") #change the current working directory to your desired diretory
os.chdir(path) #the folder to store your analysis results

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz") #download the HM27 assembly
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz") #HM46 assembly
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz") #HM65 assembly
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz") #HM69 assembly

#directory = "/home/pokoro/Paul_Okoro/"
os.system("gunzip GCF_000387825.2_ASM38782v2_genomic.fna.gz") #unzip HM27
os.system("gunzip GCF_000387845.2_ASM38784v2_genomic.fna.gz") #unzip HM46
os.system("gunzip GCF_000387785.2_ASM38778v2_genomic.fna.gz") #unzip HM65
os.system("gunzip GCF_000387865.2_ASM38786v2_genomic.fna.gz") #unzip HM69


#for file in os.listdir()
#loop through all the assembly files and count how many contigs
hm27 = []
for record in SeqIO.parse("GCF_000387825.2_ASM38782v2_genomic.fna", "fasta"):
	hm27.append(str(record.id))
open("UPEC.log", "a").write("\nThere are " + str(len(hm27)) + " contigs in the assembly HM27."+"\n")

hm46 = []
for record in SeqIO.parse("GCF_000387845.2_ASM38784v2_genomic.fna", "fasta"):
	hm46.append(str(record.id))
open("UPEC.log", "a").write("There are " + str(len(hm46)) + " contigs in the assembly HM46."+"\n")

hm65 = []
for record in SeqIO.parse("GCF_000387785.2_ASM38778v2_genomic.fna", "fasta"):
        hm46.append(str(record.id))
open("UPEC.log", "a").write("There are " + str(len(hm65)) + " contigs in the assembly HM65."+"\n")

hm69 = []
for record in SeqIO.parse("GCF_000387865.2_ASM38786v2_genomic.fna", "fasta"):
        hm69.append(str(record.id))
open("UPEC.log", "a").write("There are " + str(len(hm69)) + " contigs in the assembly HM69."+"\n")


#Loop through all the assembly files, collapse all contigs to one super contig and count length
hm27 = []
for record in SeqIO.parse("GCF_000387825.2_ASM38782v2_genomic.fna", "fasta"):
        hm27.append(str(record.seq))
open("UPEC.log", "a").write("\nThere are " + str(len(''.join(hm27))) + " bp in the assembly HM27."+"\n")

hm46 = []
for record in SeqIO.parse("GCF_000387845.2_ASM38784v2_genomic.fna", "fasta"):
        hm46.append(str(record.seq))
open("UPEC.log", "a").write("There are " + str(len(''.join(hm46))) + " bp in the assembly HM46."+"\n")

hm65 = []
for record in SeqIO.parse("GCF_000387785.2_ASM38778v2_genomic.fna", "fasta"):
        hm46.append(str(record.seq))
open("UPEC.log", "a").write("There are " + str(len(''.join(hm65))) + " bp in the assembly HM65."+"\n")

hm69 = []
for record in SeqIO.parse("GCF_000387865.2_ASM38786v2_genomic.fna", "fasta"):
        hm69.append(str(record.seq))
open("UPEC.log", "a").write("There are " + str(len(''.join(hm69))) + " bp in the assembly HM69."+"\n")

#Use prokka to annotate the assembly using the Escherichia database
os.system("prokka  --outdir HM27 --prefix HM27 --genus Escherichia --usegenus GCF_000387825.2_ASM38782v2_genomic.fna")
open("UPEC.log", "a").write("\nprokka  --outdir HM27 --prefix HM27 --genus Escherichia --usegenus GCF_000387825.2_ASM38782v2_genomic.fna"+"\n")

os.system("prokka  --outdir HM46 --prefix HM46 --genus Escherichia --usegenus GCF_000387845.2_ASM38784v2_genomic.fna")
open("UPEC.log", "a").write("prokka  --outdir HM46 --prefix HM46 --genus Escherichia --usegenus GCF_000387845.2_ASM38784v2_genomic.fna"+"\n")

os.system("prokka  --outdir HM65 --prefix HM65 --genus Escherichia --usegenus GCF_000387785.2_ASM38778v2_genomic.fna")
open("UPEC.log", "a").write("prokka  --outdir HM65 --prefix HM65 --genus Escherichia --usegenus GCF_000387785.2_ASM38778v2_genomic.fna"+"\n")

os.system("prokka  --outdir HM69 --prefix HM69 --genus Escherichia --usegenus GCF_000387865.2_ASM38786v2_genomic.fna")
open("UPEC.log", "a").write("prokka  --outdir HM69 --prefix HM69 --genus Escherichia --usegenus GCF_000387865.2_ASM38786v2_genomic.fna"+"\n")

#Write the results of the assembly to the log file
open("UPEC.log","a").write("\nAnnotation Results for HM27 \n" + (open("HM27/HM27.txt").read())+"\n")
open("UPEC.log","a").write("Annotation Results for HM46 \n" + (open("HM46/HM46.txt").read())+"\n")
open("UPEC.log","a").write("Annotation Results for HM65 \n" + (open("HM65/HM65.txt").read())+"\n")
open("UPEC.log","a").write("Annotation Results for HM69 \n" + (open("HM69/HM69.txt").read())+"\n")

#Compare the predicted CDS and tRNA from the assembly with the E. coli K-12 (NC_000913
k12_cds = 4140
k12_trna = 89

trna27 = 0
cds27 = 0
with open("HM27/HM27.txt") as hm27:
	for line in hm27:
		if line.startswith("tR"):
			trna27 = int(line.strip("\n").split(" ")[1])
		if line.startswith("CD"):
			cds27 = int(line.strip("\n").split(" ")[1])
open("UPEC.log", "a").write("\nProkka found " + str(int(cds27-k12_cds)) + " additional CDS and " + str(int(k12_trna-trna27)) + " less tRNA than the RefSeq in assembly HM27\n")

trna46 = 0
cds46 = 0
with open("HM46/HM46.txt") as hm46:
        for line in hm46:
                if line.startswith("tR"):
                        trna46 = int(line.strip("\n").split(" ")[1])
                if line.startswith("CD"):
                        cds46 = int(line.strip("\n").split(" ")[1])
open("UPEC.log", "a").write("Prokka found " + str(int(cds46-k12_cds)) + " additional CDS and " + str(int(k12_trna-trna46)) + " less tRNA than the RefSeq in assembly HM46\n")

trna65 = 0
cds65 = 0
with open("HM65/HM65.txt") as hm65:
        for line in hm65:
                if line.startswith("tR"):
                        trna65 = int(line.strip("\n").split(" ")[1])
                if line.startswith("CD"):
                        cds65 = int(line.strip("\n").split(" ")[1])
open("UPEC.log", "a").write("Prokka found " + str(int(cds65-k12_cds)) + " additional CDS and " + str(int(k12_trna-trna65)) + " less tRNA than the RefSeq in assembly HM65\n")

trna69 = 0
cds69 = 0
with open("HM69/HM69.txt") as hm69:
        for line in hm69:
                if line.startswith("tR"):
                        trna69 = int(line.strip("\n").split(" ")[1])
                if line.startswith("CD"):
                        cds69 = int(line.strip("\n").split(" ")[1])
open("UPEC.log", "a").write("Prokka found " + str(int(cds69-k12_cds)) + " additional CDS and " + str(int(k12_trna-trna69)) + " less tRNA than the RefSeq in assembly HM69\n")

#Use TopHat & Cufflinks to map the reads of a specific strain to the genome of the strain and quantify their expression, respectively
#get the trnscriptomes of the 4 E.Coli strains
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278956/SRR1278956.sra")
os.system("fastq-dump -I --split-files SRR1278956.sra") #HM27

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278960/SRR1278960.sra")
os.system("fastq-dump -I --split-files SRR1278960.sra") #HM46

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283106/SRR1283106.sra")
os.system("fastq-dump -I --split-files SRR1283106.sra") #HM65

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra")
os.system("fastq-dump -I --split-files SRR1278963.sra") #HM69

#Download the refseq assembled genome for E. coli K-12 (NC_000913) 
#os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna")

#os.system("mv NC_000913.fna EColi_k12.fa")


#Build genome index with the E.Coli k-12 (NC_000913)
#os.system("bowtie2-build NC_000913.fna EColi_K12")
#os.system("mv NC_000913.fna EColi_K12.fa")

#Build index with the HM ecoli genome
os.system("bowtie2-build HM27/HM27.fna HM27")
os.system("bowtie2-build HM46/HM46.fna HM46")
os.system("bowtie2-build HM65/HM65.fna HM65")
os.system("bowtie2-build HM69/HM69.fna HM69")

#Copy each HM genome to the same directory with the genome index
os.system("cp HM27/HM27.fna /home/pokoro/Paul_Okoro")
os.system("cp HM46/HM46.fna /home/pokoro/Paul_Okoro")
os.system("cp HM65/HM65.fna /home/pokoro/Paul_Okoro")
os.system("cp HM69/HM69.fna /home/pokoro/Paul_Okoro")

#Change the HM genome extension to .fa
os.system("mv HM27.fna HM27.fa")
os.system("mv HM46.fna HM46.fa")
os.system("mv HM65.fna HM65.fa")
os.system("mv HM69.fna HM69.fa")

#Map reads with TopHat

os.system("tophat2 -p 2 -o HM27_tophat HM27 SRR1278956_1.fastq SRR1278956_2.fastq")
os.system("tophat2 -p 2 -o HM46_tophat HM46 SRR1278960_1.fastq SRR1278960_2.fastq")
os.system("tophat2 -p 2 -o HM65_tophat HM65 SRR1283106_1.fastq SRR1283106_2.fastq")
os.system("tophat2 -p 2 -o HM69_tophat HM69 SRR1278963_1.fastq SRR1278963_2.fastq")

#Run the cufflinks using the GFF from prokka 
os.system("cufflinks -p 2 -G HM27/HM27.gff -o HM27_cufflinks_GFF HM27_tophat/accepted_hits.bam")
os.system("cufflinks -p 2 -G HM46/HM46.gff -o HM46_cufflinks_GFF HM46_tophat/accepted_hits.bam")
os.system("cufflinks -p 2 -G HM65/HM65.gff -o HM65_cufflinks_GFF HM65_tophat/accepted_hits.bam")
os.system("cufflinks -p 2 -G HM69/HM69.gff -o HM69_cufflinks_GFF HM69_tophat/accepted_hits.bam")

