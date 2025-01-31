#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:19:23 2024
@author: Frederic Labbe
This script scans the putative attachment sites (attL and attR) of a fasta sequence.
usage: python ScanAtt.py --inputfile '/path/to/file.fasta' --tRNA  --start 50000 --end 51000 --startRNA 12000
"""

import os.path
import re
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils
import argparse
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input fasta file')
parser.add_argument('-t', "--tRNA", required = True, help = 'Path to the tRNA fasta file')
parser.add_argument('-s', "--start", type = int, required = True, help = 'Start coordinate of the prophage integrase')
parser.add_argument('-e', "--end", type = int, required = True, help = 'End coordinate of the prophage integrase')
parser.add_argument('-r', "--startRNA", type = int, required = True, help = 'Start coordinate of the prophage tRNA')
args = parser.parse_args()

def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break

def ScanAtt(inputfile, tRNA, integrase_start, integrase_end, tRNA_start):
    windstep = 1
    if os.path.exists(tRNA):
        os.chdir(os.path.dirname(os.path.abspath(tRNA)))
        tRNA_fasta_sequences = SeqIO.parse(open(tRNA),'fasta')
        for tRNA_fasta in tRNA_fasta_sequences:
            tRNA_header = tRNA_fasta.id
            tRNA_sequence = tRNA_fasta.seq
        if os.path.exists(inputfile):
            ref_fasta_sequences = SeqIO.parse(open(inputfile),'fasta')
            for ref_fasta in ref_fasta_sequences:
                ref_sequence = ref_fasta.seq
            out = pd.DataFrame(columns=['Attachment_site', 'Pattern','Start','Stop','Length'])
            if tRNA_start < integrase_end:                
                for windsize in range(10, 20):
                    for subseq in chunks(tRNA_sequence, windsize, windstep):
                        if subseq in ref_sequence[integrase_end:]:
                            start = 1 + integrase_end + SeqUtils.nt_search(str(ref_sequence[integrase_end:]), subseq)[1]
                            stop = integrase_end + SeqUtils.nt_search(str(ref_sequence[integrase_end:]), subseq)[1] + windsize
                            print("Pattern:", subseq, "; Start:", start, "; Stop:", stop, "; Length:", windsize)
                            out = out.append({'Attachment_site':"attR", 'Pattern':str(subseq), 'Start':start, 'Stop':stop, 'Length':windsize}, ignore_index = True)
            else:
                for windsize in range(10, 20):
                    for subseq in chunks(tRNA_sequence, windsize, windstep):
                        if subseq in ref_sequence[:integrase_start]:
                            start = 1 + SeqUtils.nt_search(str(ref_sequence), subseq)[1]
                            stop = SeqUtils.nt_search(str(ref_sequence), subseq)[1] + windsize
                            print("Pattern:", subseq, "; Start:", start, "; Stop:", stop, "; Length:", windsize)
                            out = out.append({'Attachment_site':"attR", 'Pattern':str(subseq), 'Start':start, 'Stop':stop, 'Length':windsize}, ignore_index = True)
            out = out[out["Length"] == max(out["Length"])]
            out = out.append({'Attachment_site':"attL", 'Pattern':str(tRNA_sequence), 'Start':re.split(':|-', tRNA_header)[1], 'Stop':re.split(':|-', tRNA_header)[2], 'Length':len(tRNA_sequence)}, ignore_index = True)
            out.to_csv(tRNA.removesuffix('.fasta') + '_att' + ".csv", sep = ',', index = False)
        else:
           sys.exit('Error: provide a valid path to the input fasta file')
    else:
       sys.exit('Error: provide a valid path to the tRNA fasta file')

if __name__ == '__main__':
    ScanAtt(args.inputfile, args.tRNA, args.start, args.end, args.startRNA)
