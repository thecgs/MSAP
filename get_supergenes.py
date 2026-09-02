#!/usr/bin/env python
# coding: utf-8

import os
import argparse
from Bio import SeqIO

def tidy_name(file):
    file = os.path.basename(file).replace('.', '_').replace('-', '_')
    file = file.replace("_mafft_nucl_trimal_aln", "")
    file = file.replace("_mafft_codon_trimal_aln", "")
    file = file.replace("_mafft_prot_trimal_aln", "")
    file = file.replace("_mafft_nucl_aln", "")
    file = file.replace("_mafft_codon_aln", "")
    file = file.replace("_mafft_prot_aln", "")

    file = file.replace("_muscle_nucl_trimal_aln", "")
    file = file.replace("_muscle_codon_trimal_aln", "")
    file = file.replace("_muscle_prot_trimal_aln", "")
    file = file.replace("_muscle_nucl_aln", "")
    file = file.replace("_muscle_codon_aln", "")
    file = file.replace("_muscle_prot_aln", "")

    file = file.replace("_clustalw2_nucl_trimal_aln", "")
    file = file.replace("_clustalw2_codon_trimal_aln", "")
    file = file.replace("_clustalw2_prot_trimal_aln", "")
    file = file.replace("_clustalw2_nucl_aln", "")
    file = file.replace("_clustalw2_codon_aln", "")
    file = file.replace("_clustalw2_prot_aln", "")

    file = file.replace("_prank_nucl_trimal_aln", "")
    file = file.replace("_prank_codon_trimal_aln", "")
    file = file.replace("_prank_prot_trimal_aln", "")
    file = file.replace("_prank_nucl_aln", "")
    file = file.replace("_prank_codon_aln", "")
    file = file.replace("_prank_prot_aln", "")

    file = file.replace("_mafft_codon_trimal_fasta", "")
    file = file.replace("_muscle_codon_trimal_fasta", "")
    file = file.replace("_clustalw2_codon_trimal_fasta", "")
    file = file.replace("_prank_codon_trimal_fasta", "")

    return file


def get_supergenes(infiles, prefix):
    infiles = [file for file in infiles if os.path.getsize(file) != 0]
    supergenes = {}

    partition = open(prefix+"_partition_finder.cfg", "w")
    part_iqtree = open(prefix+'_part_iqtree.txt', 'w')
    #part_raxml = open('part_raxml.txt', 'w')

    print(f"""## ALIGNMENT FILE ##
alignment = {prefix+"_supergenes.phy"};

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = mrbayes;

## MODEL SELECCTION: AIC | AICc | BIC #
model_selection = BIC;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]    
""", file=partition)
    
    print("#nexus\nbegin sets;", file=part_iqtree)
    
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
        
        print(f"charset {tidy_name(file)} = {start_pos}-{end_pos};", file=partition)
        #print(f"DNA, {tidy_name(file)} = {start_pos}-{end_pos}", file=part_raxml)
        print(f"    charset {tidy_name(file)} = {start_pos}-{end_pos};", file=part_iqtree)
    print(f"    charpartition my_genes = :{', :'.join([tidy_name(file) for file in infiles])};\nend;", file=part_iqtree)
    print(f"#partition my_genes = {len(infiles)} : {', '.join([tidy_name(file) for file in infiles])};\n#end;", file=partition)      
    print("""
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;    
""", file=partition)

    partition.close()
    part_iqtree.close()
    #part_raxml.close()
    out = open(prefix+"_supergenes.fasta", 'w')
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
    required.add_argument('-p', '--prefix', metavar='str', required=True, 
                          help=f'A prefix output file.')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    get_supergenes(infiles=args.input, prefix=args.prefix)
    
