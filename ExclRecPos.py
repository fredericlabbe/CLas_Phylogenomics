#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 08:32:12 2023
@author: Frederic Labbe
This script creates a file listing each SNP impacted by recombination.
It uses the "importation_status.txt" files generated by ClonalFrameML, which gives the genomic regions impacted by recombination.
usage: python ExclRecPos.py --inputfile '/path/to/file.txt' --chrom "chrom1"
"""

import os.path
import re
import pandas as pd
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to input file')
parser.add_argument('-c', "--chrom", required = True, help = 'Chromosome name')
args = parser.parse_args()

def ExclRecPos(inputfile, chrom):
    if os.path.exists(inputfile):
        outputfile = inputfile.split("/")[-1].split(".txt")[0] + "_exclpos.txt"
        excl = open(inputfile, 'r')
        header = excl.readline()
        coord = list()  
        for line in excl:
            line.split()
            line = re.split(r'\t+', line)
            beg = line[1]
            end = line[2]  
            for number in range(int(beg), int(end), 1):
                coord.append(number)
        coord = list(dict.fromkeys(coord))
        df = pd.DataFrame (coord, columns = ['pos'])   
        df['chrom'] = chrom
        df = df[['chrom', 'pos']]
        df.to_csv(outputfile, sep = '\t', index = False, header = False)
    else:
        sys.exit('Error: provide a valid path to the inputfile') 

if __name__ == '__main__':
     ExclRecPos(args.inputfile, args.chrom)