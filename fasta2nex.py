#!/usr/bin/env python
# coding: utf-8

import argparse
from Bio import AlignIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="fasta format convert nex format.", add_help=False, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Date:2025/12/11 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', required=True,
                          help='A input file of fasta format.')
    required.add_argument('-o', '--output', metavar='str', required=True,
                          help='A output file of nex format.')
    optional.add_argument('-st', '--seqtype', metavar='str', default='DNA', choices=["DNA", "RNA", "protein"],
                          help=f'Sequence type. such as DNA, RNA or protein. default=DNA')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v2.00',
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    
    handle = open(args.output, 'w')
    AlignIO.convert(args.input, 'fasta', handle, 'nexus', args.seqtype)
    handle.close()
