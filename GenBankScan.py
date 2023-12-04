#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 09:55:07 2022
@author: Frederic Labbe
This script downloads genbank assemblies for a given search term (e.g. a species name).
Note: https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
usage: python GenBankScan.py --term 'Candidatus Liberibacter asiaticus' --format 'GenBank' --email 'other@example.com'
"""

import os.path
import os, shutil
from Bio import Entrez
import urllib.request
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--term", required = True, help = 'Search term, usually organism name')
parser.add_argument('-f', "--format", required = True, help = 'Accession format (RefSeq or GenBank)')
parser.add_argument('-e', "--email", required = True, help = 'Provide your own email address here')
args = parser.parse_args()

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db = "assembly", id = id, report = "full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term, format, email):
    Entrez.email = email
    handle = Entrez.esearch(db = "assembly", term = term, retmax = '200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    if not os.path.exists('old_assemblies'):
        os.mkdir('old_assemblies')
    for id in ids:
        summary = get_assembly_summary(id)
        if format == 'RefSeq':
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        elif format == 'GenBank':
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        else:
            sys.exit('Error: provide a valid accession format (RefSeq or GenBank)')
        if url == '':
            continue
        label = os.path.basename(url)        
        
        # Check file exist
        os.path.join(url, label + '_genomic.fna.gz')
        if os.path.exists(os.path.join(label + '.fna.gz')):
            print(f'The {label} assembly already exists and is readable.')
        elif os.path.exists(os.path.join(os.getcwd() + '/old_assemblies/' + label + '.fna.gz')):
            print(f'The {label} old assembly already exists and is readable.')
        else:
            #get the fasta link - change this to get other formats
            print(f'The {label} assembly does not exist and will be downloaded:')
            link = os.path.join(url, label + '_genomic.fna.gz')
            print(link)
            links.append(link)
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
            if label[-1:].isnumeric():
                version = int(label[-1:])
                newer_version = version + 1
                assembly = label[:14]
                name = label[15:-1]
                newer_assembly = assembly + str(newer_version)
                newer_name = name + str(newer_version)
                newer_label = newer_assembly + newer_name
                if os.path.exists(os.path.join(newer_label + '.fna.gz')):                            
                    print(f'A more recent version exist, the {label} assembly will be moved to the "old_assemblies" directory.')
                    shutil.move(label + '.fna.gz', 'old_assemblies/' + label + '.fna.gz')
                if version > 1:                    
                    previous_version = version - 1
                    previous_assembly = assembly + str(previous_version)
                    previous_name = name + str(previous_version)
                    previous_label = previous_assembly + previous_name
                    if os.path.exists(os.path.join(previous_label + '.fna.gz')):
                        print(f'A more recent version exist, the {previous_label} assembly will be moved to the "old_assemblies" directory.')
                        shutil.move(previous_label + '.fna.gz', 'old_assemblies/' + previous_label + '.fna.gz')
    return links

if __name__ == '__main__':
    get_assemblies(args.term, args.format, args.email)
