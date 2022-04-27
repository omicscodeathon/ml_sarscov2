#!/bin/bash
#Script to perform alignment using mafft
#COPYRIGHT 2022
#Created by Mike Mwanga
#mikemwanga6@gmail.com

#SBATCH -J Mapping
#SBATCH -p batch
#SBATCH -n 5
#SBATCH --mem-per-cpu 8000
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --mail-type=ALL

ali	
module load mafft/7.475 

mafft --auto --reorder --anysymbol sequence.fasta > sequences.aln.fas


rm sequence.fasta
mv sequences.aln.fas sequence.fasta

sbatch /home/mmwanga/ML/Blast/scripts/blast_run.sh