#!/bin/bash
#Script to extract SARS-Cov-2 genome using blast.
#COPYRIGHT 2022
#Created by Mike Mwanga
#mikemwanga6@gmail.com


#motivation from this link
#https://www.biostars.org/p/433926/

while read id start stop; do
	blastdbcmd -db sequences.fasta -entry $id -range $start-$stop
done < hits.txt
