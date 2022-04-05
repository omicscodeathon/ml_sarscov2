
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

#create a blast database
makeblastdb -in sequences.fasta -dbtype nucl -parse_seqids

#extract data for mapping spike protein
blastn -query /home/mmwanga/Blast/wuhan_spike_reference.fasta \
		-db sequences.fasta \
		-out results.out \
		-outfmt 6 \
		-max_target_seqs ${n}

#extract the sequence id, and start and stop positios
awk -F "\t" '{OFS="\t"}{print $2,$9,$10}' results.out > hits.txt


#extract a fasta file containing spike protein sequence.

#blastdbcmd -db sequences.fasta -entry MZ380280.1 -range 21563-25384

sh /home/mmwanga/Blast/scripts/read.sh > spike_protein.fas