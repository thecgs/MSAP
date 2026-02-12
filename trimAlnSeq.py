#!/usr/bin/env python
# coding: utf-8

import sys
import gzip
import argparse
from Bio import SeqIO
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Trim fasta format alignment sequence.", add_help=False, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Date:2025/09/29 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', required=True,
                          help='A input file of fasta format.')
    optional.add_argument('-o', '--output', metavar='str', default='-',
                          help='A output file of fasta format. default=None.')
    optional.add_argument('-st', '--seqtype', metavar='str', default='codon', choices=['codon', 'prot', 'nucl'],
                          help=f'Sequence type. such as nucl, prot or codon. default=codon')
    optional.add_argument('-kg', '--keepgap', action='store_true',
                          help=f'Keep the column containing gap.')
    optional.add_argument('-kn', '--keepN', action='store_true',
                          help=f'Keep the column containing the codon of N character.')
    optional.add_argument('-kx', '--keepX', action='store_true',
                          help=f'The column containing X characters is reserved. For the protein sequence.')
    optional.add_argument('-ra', '--removeall', action='store_true',
                          help=f'Remove all same columns.')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v2.00',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    
    infile = args.input
    outfile = args.output
    kg = args.keepgap
    kn = args.keepN
    kx = args.keepX
    ra = args.removeall
    st = args.seqtype
    
    indexs = []
    aln = defaultdict(list)
    if infile.endswith('.gz'):
        handle = gzip.open(infile, 'rt')
    elif infile == "-":
        handle = sys.stdin
    else:
        handle = open(infile, 'r')
    
    for record in SeqIO.parse(handle, 'fasta'):
        #print("\n"+record.id, 'length:', len(record.seq), file=sys.stderr)
        if st == 'codon':
            for n, i in enumerate(range(0, len(record.seq), 3)):
                #print(n, file=sys.stderr, end='\r')
                codon = str(record.seq.upper())[i:i+3]
                aln[record.id].append(codon)
                
                if (kn == False) or (kg == False):
                    if ('N' in codon) or ('-' in codon):
                        indexs.append(n)
                        
                if (kn == False) or (kg == True):
                    if ('N' in codon):
                        indexs.append(n)
                        
                if (kn == True) or (kg == False):
                    if ('-' in codon):
                        indexs.append(n)
                        
        elif st == 'prot':
            for n, a in enumerate(record.seq.upper()):
                #print(n, file=sys.stderr, end='\r')
                aln[record.id].append(a)
                if (kx == False) or (kg == False):
                    if ('X' == a) or ('-' == a):
                        indexs.append(n)
                        
                if (kx == False) or (kg == True):
                    if ('X' == a):
                        indexs.append(n)
                        
                if (kx == True) or (kg == False):
                    if ('-' == a):
                        indexs.append(n)
                        
        elif st == 'nucl':
            for n, a in enumerate(record.seq.upper()):
                #print(n, file=sys.stderr, end='\r')
                aln[record.id].append(a)
                if (kn == False) or (kg == False):
                    if ('N' == a) or ('-' == a):
                        indexs.append(n)
                        
                if (kn == False) or (kg == True):
                    if ('N' == a):
                        indexs.append(n)
                        
                if (kn == True) or (kg == False):
                    if ('-' == a):
                        indexs.append(n)
    handle.close()
    
    if ra == True:
        ref = list(aln.values())[0]
        for n, i in enumerate(ref):
            _status = False
            for s in list(aln.keys())[1:]:
                if aln[s][n] != ref[n]:
                    _status = True
            if _status == False:
                indexs.append(n)
                
    indexs = sorted(set(indexs))
    
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
