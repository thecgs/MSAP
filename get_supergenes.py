#!/usr/bin/env python
# coding: utf-8

import os
import argparse
from Bio import SeqIO

def get_supergenes(infiles, outfile):
    infiles = [file for file in infiles if os.path.getsize(file) != 0]
    supergenes = {}
    partition = open("partition_finder.cfg", "w")
    
    print(f"""## ALIGNMENT FILE ##
alignment = {outfile};

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = mrbayes;

## MODEL SELECCTION: AIC | AICc | BIC #
model_selection = BIC;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]    
""", file=partition)
    
    #print("#nexus\nbegin sets;", file=partition)
    
    for n, file in enumerate(infiles):
        for record in SeqIO.parse(file, 'fasta'):
            if record.id not in supergenes:
                supergenes[record.id] = record.seq
            else:
                supergenes[record.id] += record.seq
                
        if n==0:
            start_pos = 1
            end_pos = len(record.seq)
        else:
            start_pos = end_pos + 1
            end_pos = end_pos + len(record.seq)
        
        print(f"charset {os.path.basename(file).replace('.', '_')} = {start_pos}-{end_pos};", file=partition)
    
    print(f"#partition my_genes = {len(infiles)} : {', '.join([os.path.basename(file).replace('.', '_') for file in infiles])};\n#end;", file=partition)      
    print("""
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;    
""", file=partition)
    partition.close()
    out = open(outfile, 'w')
    for species in supergenes:
        print(f'>{species}\n{supergenes[species]}', file=out)
    out.close()
    return None
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get supergenes alignment file from mulitple fasta format alignment file.", 
                                     add_help=False, epilog='Date:2024/12/23 Author:Guisen Chen Email:thecgs001@foxmail.com', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='aln.fasta', 
                          help='A list of fasta format alignment files.', nargs='*', required=True)
    required.add_argument('-o', '--output', metavar='str', required=True, 
                          help=f'A output file.')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    get_supergenes(infiles=args.input, outfile=args.output)
    
