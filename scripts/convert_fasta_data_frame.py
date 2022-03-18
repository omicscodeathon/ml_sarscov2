#!/usr/bin/py
#This is a python script to convert sequence to data.frame format
#Created by Mike Mwanga
#COPYRIGHTS 2022

import sys
import re
import os.path
import fileinput


if len(sys.argv) != 3:
    print("USAGE ERROR: convert_fasta_to_data_frame.py input_file outputfile")
    sys.exit()
inputfile = sys.argv[1]
outputfile = sys.argv[2]


lines = []
with open (inputfile, 'r') as fasta_file:
    
    for line in fasta_file:
        if line.startswith(">"):
            lines.append(line.replace('\w', '').replace(">", "").replace("\n","\t"))
        next
        if not line.startswith(">"):
            lines.append(line)
with open(outputfile, 'w') as output:
    for line in lines:
        output.write(line)
