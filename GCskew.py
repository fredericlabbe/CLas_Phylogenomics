#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:31:54 2024
@author: Frederic Labbe
This script calculates the incremental and the cumulative GC skew for a given FASTA sequence according to this formula: Skew = (G - C) / (G + C)
usage: python GCskew.py --sequence '/path/to/sequence.fasta' --window 100
"""

import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC_skew
import argparse
import sys
pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser()
parser.add_argument('-s', "--sequence", required = True, help = 'Path to the FASTA sequence')
parser.add_argument('-w', "--window", type = int, required = True, help = 'Window size')
args = parser.parse_args()

def GCskew(sequence, window):
    if os.path.exists(sequence):
        phage = SeqIO.read(open(sequence), 'fasta').seq
        if len(phage) > window and 0 < window:
            gcskew = pd.DataFrame(GC_skew(phage, window = window))
            cumulative_gcskew = pd.DataFrame(np.cumsum(gcskew))
            coord = pd.DataFrame(list(range(0, len(phage), window)))
            out = pd.concat([coord, gcskew, cumulative_gcskew], axis = 1)
            out.columns =['Coordinates', 'GCskew', 'Cumulative_GCskew']
            out.to_csv(sequence.removesuffix('.fasta') + '_GCskew_Wind' + str(window) + ".csv", sep = ',', index = False)
        else:
            sys.exit('Error: provide a window size shorter than the FASTA sequence')
    else:
        sys.exit('Error: provide a valid path to the FASTA sequence')

if __name__ == '__main__':
     GCskew(args.sequence, args.window)
