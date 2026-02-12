#!/usr/bin/env python
# coding: utf-8

import os, sys
import textwrap
import argparse
from Bio import SeqIO, Seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"""Covert AA to Codon.""", 
    add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='Date:2025/01/06 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-c', '--cds', metavar='str', required=True,
                          help='A CDS file of fasta format.')
    required.add_argument('-p', '--protaln', metavar='str', required=True,
                          help='A protein alignment file of fasta format.')
    optional.add_argument('-o', '--out', metavar='str', default=None,
                          help='A CDS alignment file of fasta format.')
    optional.add_argument('-g', '--genetic_code', metavar='int', default=None,
                          type=int, help='Genetic code, check model. default=None')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v2.00',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    cdsfile=args.cds
    protalnfile = args.protaln
    outfile = args.out
    genetic_code = args.genetic_code
    
    if outfile !=None:
        out = open(outfile, 'w')
    else:
        out = sys.stdout
    
    Mapping = dict()
    for j, record in enumerate(SeqIO.parse(cdsfile, 'fasta')):
        #Mapping[record.id] = textwrap.wrap(str(record.seq.upper()), 3)
        Mapping[record.id] = textwrap.wrap(str(record.seq.upper()), 3)
    #print(Mapping)
    for j, record in enumerate(SeqIO.parse(protalnfile, 'fasta')):
        print(">"+record.id, file=out)
        s = ""
        idx = 0
        for i in record.seq.upper():
            if i != '-':
                s += Mapping[record.id][idx]
                if genetic_code!=None:
                    if str(Seq.Seq(Mapping[record.id][idx]).translate(table=genetic_code)) != i:
                        print("index:", idx, "->", str(Seq.Seq(Mapping[record.id]).translate(table=genetic_code)))
                        sys.exit()
                idx += 1
            else:
                s += '---'
        print(s, file=out)
        
    out.close()
