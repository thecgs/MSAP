## Installation

Before use, you need to install Python, and biopython.

Python3 >= 3.8

Supports four multiple sequence alignment softwares (mafft, muscle, clustalw2, prank), with mafft as the default. Therefore, the appropriate multiple sequence alignment software should be installed in the environment according to requirements before use.

```
git clone https://github.com/thecgs/MSAP.git
chmod +x *.py
pip install biopython
conda install mafft
```

Note:
1. Sequence IDs in the input FASTA format file must be unique.
2. All Python scripts require executable permissions.

## Usage

### Main Script  Usage:

```
Python MSAP.py -h
usage: MSAP.py -i str [-t int] [-s str] [-n] [-st str] [-g int] [-h] [-v]

Run multiple sequence alignment pipeline (MSAP).

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

required arguments:
  -i, --input str       A file of fasta format.

optional arguments:
  -t, --thread int      Align software thread number. default=24
  -s, --align_software str
                        Align software. such as mafft, muscle, clustalw2, prank. default=mafft
  -n, --notrim          No trim align file. The trimAlnSeq.py script used to trim codon, nucletic, or protein
                        alignment.
  -st, --seqtype str    Sequence type. such as nucl, prot or codon. default=codon
  -g, --genetic_code int
                        Genetic code, only "--model codon" take effect. default=1
  -h, --help            Show program's help message and exit.
  -v, --version         Show program's version number and exit.

Date:2026/02/12 Author:Guisen Chen Email:thecgs001@foxmail.com
```

### Format Conversion
The default format for multiple sequence alignment files is FASTA. You can use fasta2axt.py, fasta2nex.py, and fasta2phy.py to convert them to AXT, NEXUS, and PHYLIP formats respectively.
```
fasta2axt.py -i input.fasta -o output.axt
fasta2nex.py -i input.fasta -o output.nex
fasta2phy.py -i input.fasta -o output.phy
```

### Split Codon Alignment
Split the multiple sequence alignment of codons in FASTA format into positions 1, 2, and 3.
```
split_codon_seqence_alignment.py codon_aln.fasta
```

###  Construct a Super Alignment Matrix.
Merge multiple gene multiple sequence alignment files of FASTA format into a single superalignment matrix.
```
get_supergenes.py -i aln1.fasta aln2.fasta aln3.fasta -o supergene.aln.fasta
```

###  Four-Fold Degenerated Sites
Batch calculations of 4dtv (transversion rate on 4-fold degenerate sites) are performed using HKY substitution models.
```
calulate_4dtv_and_correction.py -i aln1.fasta aln2.fasta aln3.fasta -o 4dtv.tsv -g 1
```

Extract four-fold degenerate sites from a multiple sequence alignment file in FASTA format.
```
extract_4-fold_degenerated_sites.py -i aln.fasta  -o 4-fold_degenerated_sites.aln.fasta -g 1
```

### Other Script:
AA2Codon.py script is  an alternative to [pal2nal](https://github.com/liaochenlanruo/PAL2NAL)


trimAlnSeq.py script is  an alternative to [trimal](https://github.com/inab/trimal)
