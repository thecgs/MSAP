#!/usr/bin/env python
# coding: utf-8

import os
import sys
import glob
#import gzip
import argparse
import subprocess
from Bio import SeqIO, Seq
from Bio.Data import CodonTable

def remove_stop_codon(infile, outfile, genetic_code):
    out = open(outfile, 'w')
    for record in SeqIO.parse(infile, 'fasta'):
        if len(record.seq)%3 == 1:
            seqence = record.seq[:-1].upper()
        elif len(record.seq)%3 == 2:
            seqence = record.seq[:-2].upper()
        else:
            if record.seq[-3:].translate(table=genetic_code, stop_symbol='*') == '*':
                seqence = record.seq[:-3].upper()
            else:
                seqence = record.seq.upper()
        
        fix_seqence = ""   #replace inframe stop codon with NNN
        for i in range(0, len(seqence), 3):
            codon = seqence[i:i+3]
            if codon in CodonTable.unambiguous_dna_by_id[genetic_code].stop_codons:
                codon = "NNN"
            fix_seqence += codon
        print(f'>{record.id}\n{fix_seqence}', file=out)
    out.close()
    return None

def translate_seq(infile, outfile, genetic_code):
    out = open(outfile, 'w')
    for record in SeqIO.parse(infile, 'fasta'):
        if record.seq[:3] in Seq.CodonTable.unambiguous_dna_by_id[genetic_code].start_codons:
            protein_sequence = list(record.seq.translate(table=genetic_code, cds=False))
            protein_sequence[0] = "M"
            protein_sequence = ''.join(protein_sequence)
            print(f">{record.id}\n{protein_sequence}", file=out)
        else:
            print(f">{record.id}\n{record.seq.translate(table=genetic_code)}", file=out)
    out.close()
    return None

def get_prefix(infile):
    return os.path.splitext(os.path.basename(infile))[0]

def detect_seq_type(infile):
    seqence = ""
    lst = ["V", "L", "I", "E", "Q", "D", "N",
           "M", "S", "F", "W", "Y", "R", "H", 
           "P","K","X"]
    for n, record in enumerate(SeqIO.parse(infile, 'fasta')):
        if n > 1:
            break
        else:
            seqence += str(record.seq.upper())
    
    status = False
    for s in seqence:
        if s in lst:
            status = False
        else:
            if status==True:
                break
    
    if status==True:
        seq_type="PROTEIN"
    else:
        seq_type="DNA"
    return seq_type

def check_dependencies(softwares):
    #print("Check dependencies...\n")
    main_path = sys.path[0]
    
    paths = os.environ['PATH'].split(':')
    software_path = {software:None for software in softwares}
    for i in glob.glob(os.path.join(main_path, '*')):
        if os.path.isdir(i):
            paths.append(i)
    for i in glob.glob(os.path.join(main_path, '*', 'bin')):
        if os.path.isdir(i):
            paths.append(i)
    for i in glob.glob(os.path.join(main_path, '*', 'script')):
        if os.path.isdir(i):
            paths.append(i)
            
    for software in softwares:
        for path in paths:
            if os.path.exists(os.path.join(path,software)):
                software_path[software] = os.path.join(path,software)
                break
    
    for software in software_path:
        if software_path[software] == None:
            print(software, "is not installed in the environment. Please install it.")
            sys.exit()
        #else:
        #    print(software, ":", software_path[software])
    return software_path
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"""Run multiple sequence alignment pipeline (MSAP).

Codon Seqence Alignment Pipeline:
step1. Remove stop codon.
step2. The biopython used to CDS seqence convert to protein seqence.
step3. Alignment software (such as mafft (v7.525), muscle (v5.2), clustalw2 (v2.1), prank (v170427)) to align protein seqence.
step4. The AA2Codon.py script used to protein alignment convert codon aligment.
step5. The trimAlnSeq.py script used to trim codon seqence.

Nucletic or Protein Seqence Alignment Pipeline:
step1. Alignment software (such as mafft (v7.525), muscle (v5.2), clustalw2 (v2.1), prank (v170427)) to align nucletic or protein seqence.
step2. The trimAlnSeq.py script used to trim nucletic or protein seqence.

used example:
            MASP.py -i CDS.fasta -g 1 -st codon
            MASP.py -i protein.fasta -s muscle -t 1 -st prot
            MASP.py -i 16S.fasta -s mafft -t 1 -n -st nucl

Translate Tables/Genetic Codes:
1: The Standard
2: The Vertebrate Mitochondrial
3: The Yeast Mitochondrial
4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma
5: The Invertebrate Mitochondrial
6: The Ciliate, Dasycladacean and Hexamita Nuclear
9: The Echinoderm and Flatworm Mitochondrial
10: The Euplotid Nuclear
11: The Bacterial, Archaeal and Plant Plastid
12: The Alternative Yeast Nuclear
13: The Ascidian Mitochondrial
14: The Alternative Flatworm Mitochondrial
15: Blepharisma Macronuclear
16: Chlorophycean Mitochondrial
21: Trematode Mitochondrial
22: Scenedesmus obliquus Mitochondrial
23: Thraustochytrium Mitochondrial
24: Pterobranchia Mitochondrial
25: Candidate Division SR1 and Gracilibacteria
26: Pachysolen tannophilus Nuclear
27: Karyorelict Nuclear
28: Condylostoma Nuclear
29: Mesodinium Nuclear
30: Peritrich Nuclear
31: Blastocrithidia Nuclear
32: Balanophoraceae Plastid
33: Cephalodiscidae Mitochondrial
Reference website: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes

Note:
clustalw2 and prank does not allow multiple sequences to use the same name.

""", add_help=False, 
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='Date:2026/02/12 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', required=True,
                          help='A file of fasta format.')
    optional.add_argument('-t', '--thread', metavar='int', default=os.cpu_count(),
                          type=int, help=f'Align software thread number. default={os.cpu_count()}')
    optional.add_argument('-s', '--align_software', metavar='str', default='mafft',
                          choices=["mafft", "muscle", "prank", "clustalw2"], help='Align software. such as mafft, muscle, clustalw2, prank. default=mafft')    
    optional.add_argument('-n', '--notrim', action='store_true',
                          help=f'No trim align file. The trimAlnSeq.py script used to trim codon, nucletic, or protein alignment.')
    optional.add_argument('-st', '--seqtype', metavar='str', default='codon', choices=['codon', 'prot', 'nucl'],
                          help=f'Sequence type. such as nucl, prot or codon. default=codon')      
    optional.add_argument('-g', '--genetic_code', metavar='int', default=1,
                          type=int, help='Genetic code, only "--model codon" take effect. default=1')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v2.00',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    infile=args.input
    thread = args.thread
    genetic_code = args.genetic_code
    align_software = args.align_software
    seqtype = args.seqtype
    notrim= args.notrim

software_path = check_dependencies(softwares=[align_software])

# run main
prefix = get_prefix(infile)
infile = os.path.realpath(infile)
if seqtype == "codon":
    infile_tmp = prefix + '.CDS.fasta'
    remove_stop_codon(infile, outfile=infile_tmp, genetic_code=genetic_code)      # remove stop codon, if there is stop codon, codeml can not work.
    infile = infile_tmp
    translate_seq(infile, outfile=prefix + '.pep.fasta', genetic_code=genetic_code)
    align_infile = os.path.realpath(f"{prefix}.pep.fasta")
else:
    align_infile = infile
    
## run align
if align_software == "mafft":
    if seqtype == "codon" or seqtype == "prot":
        cmd = f"{software_path['mafft']} --thread {thread} --quiet --auto {align_infile} > {prefix}.mafft.prot.aln"
        align_outfile = os.path.realpath(f"{prefix}.mafft.prot.aln")
    elif seqtype == "nucl":
        cmd = f"{software_path['mafft']} --thread {thread} --quiet --auto {align_infile} > {prefix}.mafft.nucl.aln"
        align_outfile = os.path.realpath(f"{prefix}.mafft.nucl.aln")
    subprocess.run(cmd, shell=True, capture_output=subprocess.PIPE)
    
    
elif align_software == "muscle":
    if seqtype == "codon" or seqtype == "prot":
        cmd = f"{software_path['muscle']} -threads {thread} -align {align_infile} -output {prefix}.muscle.prot.aln"
        align_outfile = os.path.realpath(f"{prefix}.muscle.prot.aln")
    elif seqtype == "nucl":
        cmd = f"{software_path['muscle']} -threads {thread} -align {align_infile} -output {prefix}.muscle.nucl.aln"
        align_outfile = os.path.realpath(f"{prefix}.muscle.nucl.aln")
    subprocess.run(cmd, shell=True, capture_output=subprocess.PIPE)
    
    
elif align_software == "prank":
    cmd = f"{software_path['prank']} -d={align_infile} -o={prefix}.prank.aln -f=fasta"
    subprocess.run(cmd, shell=True, capture_output=subprocess.PIPE)
    if seqtype == "codon" or seqtype == "prot":
        os.rename(f"{prefix}.prank.aln.best.fas", f"{prefix}.prank.prot.aln")
        align_outfile = os.path.realpath(f"{prefix}.prank.prot.aln")
    elif seqtype == "nucl":
        os.rename(f"{prefix}.prank.aln.best.fas", f"{prefix}.prank.nucl.aln")
        align_outfile = os.path.realpath(f"{prefix}.prank.nucl.aln")
        
        
elif align_software == "clustalw2":
    seq_type = detect_seq_type(infile)
    if seqtype == "codon" or seqtype == "prot":
        cmd = f"{software_path['clustalw2']} -INFILE={align_infile} -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE={prefix}.clustalw2.pep.aln"
        subprocess.run(cmd, shell=True, capture_output=subprocess.PIPE)
        if seqtype == "codon":
            os.remove(os.path.join(os.path.dirname(infile), f"{prefix}.pep.dnd"))
        else:
            os.remove(os.path.join(os.path.dirname(infile), f"{prefix}.dnd"))
        os.rename(f"{prefix}.clustalw2.pep.aln", f"{prefix}.clustalw2.prot.aln")
        align_outfile = os.path.realpath(f"{prefix}.clustalw2.prot.aln")
    if seqtype == "nucl":
        cmd = f"{software_path['clustalw2']} -INFILE={align_infile} -TYPE=DNA -OUTPUT=FASTA -OUTFILE={prefix}.clustalw2.nucl.aln"
        subprocess.run(cmd, shell=True, capture_output=subprocess.PIPE)
        os.remove(os.path.join(os.path.dirname(infile), f"{prefix}.dnd"))
        align_outfile = os.path.realpath(f"{prefix}.clustalw2.nucl.aln")
else:
    print("Please input a align software. [matff, muscle, prank, clustalw2]")
    sys.exit()
    
## run AA2codon.py
if seqtype == "codon":
    cmd = f"{os.path.join(sys.path[0], 'AA2Codon.py')} -c {infile} -p {prefix}.{align_software}.prot.aln -o {prefix}.{align_software}.codon.aln"
    #subprocess.run(cmd, shell=True, capture_output=True)
    subprocess.run(cmd, shell=True)
    if notrim == False:
        cmd = f"{os.path.join(sys.path[0], 'trimAlnSeq.py')} -i  {prefix}.{align_software}.codon.aln -o {prefix}.{align_software}.codon.trimal.aln -st codon"
        subprocess.run(cmd, shell=True, capture_output=True)
else:
    ## run trimAlnSeq.py
    if notrim == False:
        cmd = f"{os.path.join(sys.path[0], 'trimAlnSeq.py')} -i  {prefix}.{align_software}.{seqtype}.aln -o {prefix}.{align_software}.{seqtype}.trimal.aln -st {seqtype}"
        subprocess.run(cmd, shell=True)
        
if seqtype == "codon":
    os.remove(infile)
    

### run pal2nal
#if seqtype == "codon":
#    cmd = f"{pal2nal} {prefix}.{align_software}.prot.aln {infile} -codontable {genetic_code} -output fasta > {prefix}.{align_software}.codon.aln"
#    subprocess.run(cmd, shell=True, capture_output=True)
#    if notrim == False:
#        cmd = f"{pal2nal} {prefix}.{align_software}.prot.aln {infile} -codontable {genetic_code} -output fasta -nogap > {prefix}.{align_software}.codon.trimal.aln"
#        subprocess.run(cmd, shell=True, capture_output=True)
#else:
#    ## run trimal
#    if notrim == False:
#        cmd = f"{trimal} -in {prefix}.{align_software}.{seqtype}.aln -out {prefix}.{align_software}.{seqtype}.trimal.aln -automated1"
#        subprocess.run(cmd, shell=True)
#        
#if seqtype == "codon":
#    os.remove(infile)
