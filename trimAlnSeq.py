#!/usr/bin/env python
# coding: utf-8

import sys
import gzip
import argparse
from Bio import SeqIO
from collections import defaultdict, Counter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Trim fasta format alignment sequence.", add_help=False, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Date:2026/06/21 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', required=True,
                          help='A input file of fasta format.')
    optional.add_argument('-o', '--output', metavar='str', default='-',
                          help='A output file of fasta format. default=None')
    optional.add_argument('-st', '--seqtype', metavar='str', default='codon', choices=['codon', 'prot', 'nucl'],
                          help=f'Sequence type. such as nucl, prot or codon. default=codon')
    
    optional.add_argument('-G', '--G', default=0, type=float, metavar='float',
                          help=f'Gap maxinum ratio for per site. range 0-1. default=0')
    optional.add_argument('-N', '--N', default=0, type=float, metavar='float',
                          help=f'N maxinum ratio for per site in Nucl or Codon sequence. range 0-1. default=0')
    optional.add_argument('-X', '--X', default=0, type=float, metavar='float',
                          help=f'X maxinum ratio for per site in protein sequence. range 0-1. default=0')
    
    optional.add_argument('-ra', '--removeall', action='store_true',
                          help=f'Remove all same columns.')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v2.01',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    
    infile = args.input
    outfile = args.output
    ra = args.removeall
    st = args.seqtype
    G = args.G
    N = args.N
    X = args.X
    
    if infile.endswith('.gz'):
        handle = gzip.open(infile, 'rt')
    elif infile == "-":
        handle = sys.stdin
    else:
        handle = open(infile, 'r')
    
    ## create index
    indexs = []
    aln = defaultdict(list)
    if st!='codon':
        Index2Seq = defaultdict(list)
        for record in SeqIO.parse(handle, 'fasta'):
            for n, s in enumerate(str(record.seq.upper())):
                aln[record.id].append(s)
                Index2Seq[n].append(s)
        handle.close()
    else:
        Index2Seq = defaultdict(list)
        for record in SeqIO.parse(handle, 'fasta'):
            for n, i in enumerate(range(0, len(record.seq), 3)):
                codon = str(record.seq.upper())[i:i+3]
                aln[record.id].append(codon)
                if 'N' in codon:
                    Index2Seq[n].append("N")
                elif '-' in codon:
                    Index2Seq[n].append("-")
                else:
                    Index2Seq[n].append(codon)
        handle.close()
    
    #print(Index2Seq)
    ## get filter site index
    if st!='prot':
        seqlen = len(Index2Seq.keys())
        seqnum = len(Index2Seq[0])
        for n in range(0, seqlen):
            N_Ratio = Counter(Index2Seq[n]).get("N", 0) / seqnum
            gap_Ratio = Counter(Index2Seq[n]).get("-", 0) / seqnum
            if N < N_Ratio:
                indexs.append(n)
            if G < gap_Ratio:
                indexs.append(n)
            if (ra == True) and (len(Counter(Index2Seq[n]))) == 1:
                indexs.append(n)
    else:
        seqlen = len(Index2Seq.keys())
        seqnum = len(Index2Seq[0])
        for n in range(0, seqlen):
            X_Ratio = Counter(Index2Seq[n]).get("X", 0) / seqnum
            gap_Ratio = Counter(Index2Seq[n]).get("-", 0) / seqnum
            if X < X_Ratio:
                indexs.append(n)
            if G < gap_Ratio:
                indexs.append(n)
            if (ra == True) and (len(Counter(Index2Seq[n]))) == 1:
                indexs.append(n)
                
    indexs = sorted(set(indexs))
    #print(indexs)
    if len(indexs) == len(list(aln.values())[0]):
        sys.exit()
        
    if outfile == "-":
        out = sys.stdout
    elif outfile.endswith('.gz'):
        out = gzip.open(outfile, 'wt')
    else:
        out = open(outfile, 'w')
        
    for s in aln:
        print(">"+s, file=out)
        seqs = ""
        for n, codon in enumerate(aln[s]):
            if n not in indexs:
                seqs += codon
        print(seqs, file=out)
    out.close()
