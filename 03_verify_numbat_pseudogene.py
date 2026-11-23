{\rtf1\ansi\ansicpg1252\cocoartf2908
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/env python3\
"""\
Step 3: Audit all 45 marsupial TAS1R sequences for internal stop codons and extract\
61-bp genomic probes spanning the two nonsense mutations in Myrmecobius fasciatus TAS1R1.\
"""\
\
import os\
from Bio import SeqIO\
\
LOCI_FNA = os.path.expanduser("~/tas1r_pipeline/results/cds/tas1r_genomic_loci.fna")\
CDS_FNA = os.path.expanduser("~/tas1r_pipeline/results/cds/tas1r_marsupials_cds.fna")\
\
# Extract 61-bp genomic probes for Myrmecobius fasciatus TAS1R1 stop codons\
locus = next(r for r in SeqIO.parse(LOCI_FNA, "fasta") if r.id.startswith("Myrmecobius_fasciatus|TAS1R1"))\
seq_str = str(locus.seq)\
\
print(">Myrmecobius_TAS1R1_Stop1_Exon2_TAG_codon148")\
print(seq_str[3479:3540])\
print(">Myrmecobius_TAS1R1_Stop2_Exon4_TGA_codon467")\
print(seq_str[6791:6852])\
\
# Full stop-codon audit across all 45 marsupial CDS sequences\
for rec in SeqIO.parse(CDS_FNA, "fasta"):\
    s = str(rec.seq)\
    trim = len(s) - (len(s) % 3)\
    aa = str(rec.seq[:trim].translate(to_stop=False))\
    body_aa = aa[:-1] if aa.endswith("*") else aa\
    if "*" in body_aa:\
        stops = [i + 1 for i, c in enumerate(body_aa) if c == "*"]\
        print(f"PSEUDOGENE: \{rec.id\} -> internal stop codons at aa positions \{stops\}")\
}