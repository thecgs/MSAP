#!/usr/bin/env python
# coding: utf-8

import argparse
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import Counter, defaultdict

def get_four_fold_codons(genetic_code=1):
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]
    #print('The condon table used:', codon_table, file=sys.stderr)
    four_fold_codons = {}
    acids2codons = defaultdict(list)
    for codon in codon_table.forward_table:
        acids2codons[codon_table.forward_table[codon]].append(codon)
        
    for acid in acids2codons:
        tmp = defaultdict(list)
        for codon in acids2codons[acid]:
            tmp[codon[0:2]].append(codon)
        for i in tmp:
            if len(tmp[i]) == 4:
                for c in tmp[i]:
                    four_fold_codons.setdefault(c, acid)
    
    return four_fold_codons

def extract_four_fold_site(infile, outfile, genetic_code):
    alignment_seqences = {}
    for record in SeqIO.parse(infile, 'fasta'):
        alignment_seqences.setdefault(record.id, record.seq.upper())
    
    four_fold_codons = get_four_fold_codons(genetic_code)
    alignment_seqences_four_fold_site = defaultdict(str)

    for i in range(0, len(record.seq), 3):
        #print(i+3, end='\r')
        codons = {}
        for geneID in alignment_seqences:
            codons[geneID] = alignment_seqences[geneID][i: i+3]

        if len(Counter([c[0:2] for c in codons.values()])) == 1:
            codons_bool = [True if c in four_fold_codons else False for c in codons.values()]
            if all(codons_bool):
                for geneID in alignment_seqences:
                    alignment_seqences_four_fold_site[geneID] += codons[geneID]
                    
    out = open(outfile, 'w')
    for k in alignment_seqences_four_fold_site:
        print(f'>{k}\n{alignment_seqences_four_fold_site[k]}', file=out)
    out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""extract four-fold degenerated sites

Translate Tables/Genetic Codes:
 1: Standard
 2: Vertebrate Mitochondrial
 3: YeastMitochondrial
 4: Mold Mitochondrial, Protozoan Mitochondrial, Coelenterate Mitochondrial, Mycoplasma, Spiroplasma
 5: Invertebrate Mitochondrial
 6: Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear
 9: Echinoderm Mitochondrial, Flatworm Mitochondrial
10: Euplotid Nuclear
11: Bacterial, Archaeal, Plant Plastid
12: Alternative Yeast Nuclear
13: Ascidian Mitochondrial
14: Alternative Flatworm Mitochondrial
16: Chlorophycean Mitochondrial
21: Trematode Mitochondrial
22: Scenedesmus obliquus Mitochondrial
23: Thraustochytrium Mitochondrial
24: Rhabdopleuridae Mitochondrial
25: Candidate Division SR1, Gracilibacteria
26: Pachysolen tannophilus Nuclear
27: Karyorelict Nuclear
28: Condylostoma Nuclear
29: Mesodinium Nuclear
30: Peritrich Nuclear
31: Blastocrithidia Nuclear
33: Cephalodiscidae Mitochondrial UAA-Tyr

""", add_help=False, epilog='Date:2024/12/25 Author:Guisen Chen Email:thecgs001@foxmail.com', formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', help='A fasta format input file.', required=True)    
    required.add_argument('-o', '--output', metavar='str', help='A fasta format output file.', required=True)
    optional.add_argument('-g', '--genetic_code', metavar='int', default=1, type=int, help=f'Genetic code. default=1')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    extract_four_fold_site(infile=args.input, outfile=args.output, genetic_code=args.genetic_code)
