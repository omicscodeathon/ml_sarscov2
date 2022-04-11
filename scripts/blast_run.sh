#!/bin/bash
#Script to extract spike protein from complete SARS-Cov-2 genome using blast.
#COPYRIGHT 2022
#Created by Mike Mwanga
#mikemwanga6@gmail.com
#motivation from this link
#https://www.biostars.org/p/433926/

module load blast/2.7.1+ 

#get number of sequences
n=$(grep -c ">" sequences.fasta)

echo create database **********
#create a blast database
makeblastdb -in sequences.fasta -dbtype nucl -parse_seqids
 
echo ****Running blast******
#extract data for mapping spike protein
blastn -query /home/mmwanga/Blast/wuhan_spike_reference.fasta -db sequences.fasta -out results.out -outfmt 6 -max_hsps 1 -max_target_seqs ${n}

echo ****Extracting coordinates******

#extract the sequence id, and start and stop positios
awk -F "\t" '{OFS="\t"}{print $2,$9,$10}' results.out > hits.txt

python3 /home/mmwanga/Blast/scripts/extract.py

echo ***Getting the sequences******
sh /home/mmwanga/Blast/scripts/read.sh > spike_protein.fas

echo ****Run alignment*****

sbatch /home/mmwanga/Blast/scripts/alignment.sh

echo **Completed**


