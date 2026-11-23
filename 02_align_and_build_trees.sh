{\rtf1\ansi\ansicpg1252\cocoartf2908
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/env bash\
set -euo pipefail\
\
BASE_DIR="$HOME/tas1r_pipeline"\
ALN_DIR="$BASE_DIR/results/alignments"\
TREE_DIR="$BASE_DIR/results/trees"\
mkdir -p "$ALN_DIR" "$TREE_DIR"\
\
for gene in TAS1R1 TAS1R2 TAS1R3; do\
    # 1. Translate 18-taxon CDS (15 marsupials + 3 placental outgroups) to amino acids\
    python3 - << EOF\
from Bio import SeqIO\
from Bio.Seq import Seq\
\
cds_in = "$BASE_DIR/results/cds/$\{gene\}_18taxa.fna"\
prot_out = "$ALN_DIR/$\{gene\}_18taxa_proteins.faa"\
with open(prot_out, "w") as out:\
    for rec in SeqIO.parse(cds_in, "fasta"):\
        s = str(rec.seq).upper()\
        trim = len(s) - (len(s) % 3)\
        aa = str(Seq(s[:trim]).translate(to_stop=False)).rstrip("*")\
        out.write(f">\{rec.id\}\\n\{aa\}\\n")\
EOF\
\
    # 2. Align amino acid sequences using MAFFT (L-INS-i)\
    mafft --auto "$ALN_DIR/$\{gene\}_18taxa_proteins.faa" > "$ALN_DIR/$\{gene\}_18taxa_protein_aln.faa"\
\
    # 3. Thread nucleotide triplets back onto the protein alignment (frame-preserving)\
    python3 - << EOF\
from Bio import SeqIO\
from Bio.Seq import Seq\
\
prot_aln = "$ALN_DIR/$\{gene\}_18taxa_protein_aln.faa"\
cds_in = "$BASE_DIR/results/cds/$\{gene\}_18taxa.fna"\
cod_out = "$ALN_DIR/$\{gene\}_18taxa_codon_aln.fna"\
\
cds_dict = \{r.id: str(r.seq).upper() for r in SeqIO.parse(cds_in, "fasta")\}\
with open(cod_out, "w") as out:\
    for prec in SeqIO.parse(prot_aln, "fasta"):\
        cds = cds_dict[prec.id]\
        idx = 0\
        aligned_nt = []\
        for aa in str(prec.seq):\
            if aa == "-":\
                aligned_nt.append("---")\
            else:\
                codon = cds[idx:idx+3]\
                if len(codon) == 3 and str(Seq(codon).translate()) != aa:\
                    if idx + 4 <= len(cds) and str(Seq(cds[idx+1:idx+4]).translate()) == aa:\
                        idx += 1\
                        codon = cds[idx:idx+3]\
                    elif idx + 5 <= len(cds) and str(Seq(cds[idx+2:idx+5]).translate()) == aa:\
                        idx += 2\
                        codon = cds[idx:idx+3]\
                aligned_nt.append(codon if len(codon) == 3 else "---")\
                idx += 3\
        out.write(f">\{prec.id\}\\n\{''.join(aligned_nt)\}\\n")\
EOF\
\
    # 4. Infer rooted Maximum-Likelihood phylogeny with ModelFinder + 1,000 UFBoot2 replicates\
    iqtree -s "$ALN_DIR/$\{gene\}_18taxa_codon_aln.fna" \\\
        -m MFP -bb 1000 -nt AUTO \\\
        -o Homo_sapiens,Mus_musculus,Loxodonta_africana \\\
        -pre "$TREE_DIR/$\{gene\}_18taxa" -redo\
done\
}