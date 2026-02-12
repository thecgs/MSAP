#!/usr/bin/env python
# coding: utf-8

import os
import sys
import math
import argparse
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict

def get_four_fold_codons(genetic_code=1):
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]
    print('The condon table used:', codon_table, file=sys.stderr)
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

def get_len(seq):
    if len(seq)%3 == 0:
        l = len(seq)
    elif len(seq)%3 == 1:
        l = len(seq) - 1
    else:
        l = len(seq) - 2
    return l

def is_transversion(base1, base2):
    transversion = {'A': 'TC', 'G': 'TC', 'C': 'AG', 'T':'AG'}
    if base2 in transversion[base1]:
        res = True
    else:
        res = False
    return res

def calculate_4DTV_correction(infile, four_fold_codons, error=sys.stderr):
    print(infile)
    file_prefix = os.path.splitext(os.path.basename(infile))[0]
    
    for i, record in enumerate(SeqIO.parse(infile, 'fasta')):
        if i == 0:
            seq1 = record.seq.upper()
            seq1_name = record.id
        elif i == 1:
            seq2 = record.seq.upper()
            seq2_name = record.id
        else:
            print(f'Warning: The {infile} alignment file contains more than two sequences, only the first two records are used.', file=error)
            break
            #sys.exit()
    
    #print(infile, len(seq1), len(seq2))
    if len(seq1) != len(seq2):
        print('Error: sequence1 length is not equal to sequence2 length, please check if it is a alignment file.', file=error)
        sys.exit()
    
    seq = ''
    fourfold_sites_total_number = 0
    fourfold_sites_transversion_number = 0
    for i in range(0, get_len(seq1), 3):
        codon1 = seq1[i: i+3]
        codon2 = seq2[i: i+3]
        if (codon1 in four_fold_codons) and (codon2 in four_fold_codons) and (codon1[0:2] == codon2[0:2]):
            fourfold_sites_total_number += 1
            seq += codon1[2]
            seq += codon2[2]
            if is_transversion(codon1[2], codon2[2]):
                fourfold_sites_transversion_number += 1

    raw_4dtv = fourfold_sites_transversion_number/fourfold_sites_total_number

    A = 0.5*seq.count('A')/fourfold_sites_total_number
    C = 0.5*seq.count('C')/fourfold_sites_total_number
    G = 0.5*seq.count('G')/fourfold_sites_total_number
    T = 0.5*seq.count('T')/fourfold_sites_total_number
    Y = 0.5*(seq.count('T') + seq.count('C'))/fourfold_sites_total_number
    R = 0.5*(seq.count('A') + seq.count('G'))/fourfold_sites_total_number

    if (A != 0) and (C != 0) and (G != 0) and (T != 0) and (Y != 0) and (R != 0):
        if (1-raw_4dtv*(T*C*R/Y+A*G*Y/R)/(2*(T*C*R+A*G*Y))) == 0:
            corrected_4dtv = "NA"
        else:
            a= -1*math.log(1-raw_4dtv*(T*C*R/Y+A*G*Y/R)/(2*(T*C*R+A*G*Y)))
            if (1-raw_4dtv/(2*Y*R) > 0):
                b=-1*math.log(1-raw_4dtv/(2*Y*R))
                corrected_4dtv=2*a*(T*C/Y+A*G/R)-2*b*(T*C*R/Y+A*G*Y/R-Y*R)
            else:
                corrected_4dtv = "NA"
    else:
        corrected_4dtv = 'NA'
    return file_prefix, seq1_name, seq2_name, corrected_4dtv, raw_4dtv, fourfold_sites_total_number, fourfold_sites_transversion_number

def main(infiles, outfile, genetic_code=1):
    four_fold_codons = get_four_fold_codons(genetic_code=genetic_code)
    #print(four_fold_codons)
    out = open(outfile, 'w')
    print('Input file prefix\tSeqence1 name\tSeqence2 name\tcorrected_4dtv\traw_4dtv\tfourfold_sites_total_number\tfourfold_sites_transversion_number', file=out)
    for infile in infiles:
        file_prefix, seq1_name, seq2_name, corrected_4dtv, raw_4dtv, fourfold_sites_total_number, fourfold_sites_transversion_number = calculate_4DTV_correction(infile, four_fold_codons)
        print(file_prefix, seq1_name, seq2_name, corrected_4dtv, raw_4dtv, fourfold_sites_total_number, fourfold_sites_transversion_number, sep='\t', file=out)
    out.close()
    return None
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""4dtv (transversion rate on 4-fold degenerated sites) are calculated with HKY substitution models 

Reference: M. Hasegawa, H. Kishino, and T. Yano, J. Mol. Evol. 22 (2), 160 (1985)

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

""", add_help=False, epilog='Date:2024/12/24 Author:Guisen Chen Email:thecgs001@foxmail.com', formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', help='A fasta format input file.', required=True, nargs='*')    
    required.add_argument('-o', '--output', metavar='str', help='A tsv format output file.', required=True)
    optional.add_argument('-g', '--genetic_code', metavar='int', default=1, type=int, help=f'Genetic code. default=1')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    main(infiles=args.input, outfile=args.output, genetic_code=args.genetic_code)
