#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from Bio import SeqIO

def fasta2axt(inputfile, output=None):
    ID = []
    Seq = []
    for record in SeqIO.parse(inputfile, 'fasta'):
        ID.append(record.id)
        Seq.append(str(record.seq))
    
    if output ==None:
        axt_out = sys.stdout
    else:
        axt_out = open(output, 'w')
        
    print(">"+"|".join(ID), file=axt_out)
    print("\n".join(Seq), file=axt_out)
    axt_out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fasta format convert axt format.',  add_help=False,
                                     epilog='Date:2025/03/06 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', 
                          help='A input file of fasta format.')
    optional.add_argument('-o', '--output', metavar='str', 
                          help='A output file of axt format. (default=sys.stdout)', default=None)
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00', 
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    fasta2axt(inputfile=args.input, output=args.output)
