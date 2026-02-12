#!/usr/bin/env python
# coding: utf-8

import os
import argparse
from Bio import SeqIO

def get_prefix(infile):
    return os.path.splitext(os.path.basename(infile))[0]

def split(infile):

    out0 = open('codon1.' + get_prefix(infile) + '.fasta', 'w')
    out1 = open('codon2.' + get_prefix(infile) + '.fasta', 'w')
    out2 = open('codon3.' + get_prefix(infile) + '.fasta', 'w')

    for record in SeqIO.parse(infile, 'fasta'):
        codon1 = record.seq[::3]
        codon2 = record.seq[1::3]
        codon3 = record.seq[2::3]

        print(f">{record.id}\n{codon1}", file=out0)
        print(f">{record.id}\n{codon2}", file=out1)
        print(f">{record.id}\n{codon3}", file=out2)

    out0.close()
    out1.close()
    out2.close()
    
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='to split codon seqence alignment.',  add_help=False,
                                     epilog='date:2024/12/11 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('alnfile', metavar='codon_aln.fasta', 
                          help='A Fasta format file for codon sequence alignment')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    split(infile=args.alnfile)
