#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from Bio import SeqIO

def fasta2phy(infile, outfile):

    if outfile == None:
        out = sys.stdout
    else:
        out = open(outfile, 'w')
    
    seqences = []
    for record in SeqIO.parse(infile, 'fasta'):
        seqences.append(record)
        
    print("{} {}".format(len(seqences),len(seqences[0])), file=out)
    for record in seqences:
        print(record.id+' ', record.seq.upper(), file=out)
    out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fasta alignment convert phy alignment.', add_help=False,
                                     epilog='Date:2024/12/10 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', help='A input file of fasta format.', required=True)
    optional.add_argument('-o', '--output', metavar='str', help='A output file of phy format. default=sys.stdout', default=None)
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    fasta2phy(infile=args.input, outfile=args.output)
