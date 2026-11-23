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
HYPHY_DIR="$BASE_DIR/results/hyphy"\
mkdir -p "$HYPHY_DIR"\
\
# 1. Standardize species binomials and mask stop codons (including numbat TAS1R1 stops) as ---\
python3 - << 'EOF'\
import os\
from Bio import SeqIO\
\
NAME_MAP = \{\
    "Dasyurus_vivverinus": "Dasyurus_viverrinus",\
    "GCF_003287225.1_phaCin_HiC_genomic": "Phascolarctos_cinereus",\
    "GCF_011100635.1_mTriVul1.pri_genomic": "Trichosurus_vulpecula",\
    "GCF_016432865.1_AdamAnt_v2_genomic_Antechinus": "Antechinus_flavipes",\
    "GCF_027887165.2_mMonDom1.pri_genomic": "Monodelphis_domestica",\
    "GCF_028372415.1_mMacEug1.pri_v2_genomic": "Notamacropus_eugenii",\
    "GCF_037893015.1_bilby.v1.9.chrom.fasta_genomic": "Macrotis_lagotis",\
    "GCF_048593235.1_ASM4859323v1_genomic_Sminthopsis": "Sminthopsis_crassicaudata",\
    "GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic": "Vombatus_ursinus",\
    "GCF_902635505.1_mSarHar1.11_genomic": "Sarcophilus_harrisii",\
    "Myrmecobius_fasciatus": "Myrmecobius_fasciatus",\
    "Perameles_gunnii": "Perameles_gunnii",\
    "Thylacinus_cynocephalus": "Thylacinus_cynocephalus",\
    "kangaroo_genomic": "Macropus_giganteus",\
    "notoryctes_genomic": "Notoryctes_typhlops",\
    "Homo_sapiens": "Homo_sapiens",\
    "Mus_musculus": "Mus_musculus",\
    "Loxodonta_africana": "Loxodonta_africana",\
\}\
STOPS = \{"TAA", "TAG", "TGA"\}\
\
for gene in ("TAS1R1", "TAS1R2", "TAS1R3"):\
    in_aln = os.path.expanduser(f"~/tas1r_pipeline/results/alignments/\{gene\}_18taxa_codon_aln.fna")\
    out_aln = os.path.expanduser(f"~/tas1r_pipeline/results/hyphy/\{gene\}_18taxa_hyphy.fna")\
    cleaned = []\
    for rec in SeqIO.parse(in_aln, "fasta"):\
        sp = NAME_MAP.get(rec.id, rec.id)\
        seq = str(rec.seq).upper()\
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]\
        masked = ["---" if c in STOPS else c for c in codons]\
        cleaned.append((sp, masked))\
\
    n_codons = len(cleaned[0][1])\
    keep_cols = [i for i in range(n_codons) if any(codons[i] != "---" for _, codons in cleaned)]\
    with open(out_aln, "w") as out:\
        for sp, codons in cleaned:\
            out.write(f">\{sp\}\\n\{''.join(codons[i] for i in keep_cols)\}\\n")\
EOF\
\
# 2. Write species trees for Hypothesis 1 (Numbat foreground), Hypothesis 2 (Dietary guilds), and MEME\
cat << 'EOF' > "$HYPHY_DIR/species_tree_numbat_test.nwk"\
(((((((Dasyurus_viverrinus\{Reference\},Sarcophilus_harrisii\{Reference\})\{Reference\},Antechinus_flavipes\{Reference\})\{Reference\},Sminthopsis_crassicaudata\{Reference\})\{Reference\},Thylacinus_cynocephalus\{Reference\})\{Reference\},Myrmecobius_fasciatus\{Test\})\{Reference\},(((Phascolarctos_cinereus\{Reference\},Vombatus_ursinus\{Reference\})\{Reference\},(Trichosurus_vulpecula\{Reference\},(Notamacropus_eugenii\{Reference\},Macropus_giganteus\{Reference\})\{Reference\})\{Reference\})\{Reference\},((Macrotis_lagotis\{Reference\},Perameles_gunnii\{Reference\})\{Reference\},Notoryctes_typhlops\{Reference\})\{Reference\})\{Reference\})\{Reference\},Monodelphis_domestica\{Reference\},((Homo_sapiens,Mus_musculus),Loxodonta_africana));\
EOF\
\
cat << 'EOF' > "$HYPHY_DIR/species_tree_dietary_guilds.nwk"\
(((((((Dasyurus_viverrinus\{Test\},Sarcophilus_harrisii\{Test\})\{Test\},Antechinus_flavipes\{Test\})\{Test\},Sminthopsis_crassicaudata\{Test\})\{Test\},Thylacinus_cynocephalus\{Test\})\{Test\},Myrmecobius_fasciatus\{Test\})\{Test\},(((Phascolarctos_cinereus\{Reference\},Vombatus_ursinus\{Reference\})\{Reference\},(Trichosurus_vulpecula\{Reference\},(Notamacropus_eugenii\{Reference\},Macropus_giganteus\{Reference\})\{Reference\})\{Reference\})\{Reference\},((Macrotis_lagotis\{Reference\},Perameles_gunnii\{Reference\})\{Reference\},Notoryctes_typhlops\{Test\})\{Reference\})\{Reference\})\{Reference\},Monodelphis_domestica\{Reference\},((Homo_sapiens,Mus_musculus),Loxodonta_africana));\
EOF\
\
cat << 'EOF' > "$HYPHY_DIR/species_tree_marsupials_only.nwk"\
(((((((Dasyurus_viverrinus\{Test\},Sarcophilus_harrisii\{Test\})\{Test\},Antechinus_flavipes\{Test\})\{Test\},Sminthopsis_crassicaudata\{Test\})\{Test\},Thylacinus_cynocephalus\{Test\})\{Test\},Myrmecobius_fasciatus\{Test\})\{Test\},(((Phascolarctos_cinereus\{Test\},Vombatus_ursinus\{Test\})\{Test\},(Trichosurus_vulpecula\{Test\},(Notamacropus_eugenii\{Test\},Macropus_giganteus\{Test\})\{Test\})\{Test\})\{Test\},((Macrotis_lagotis\{Test\},Perameles_gunnii\{Test\})\{Test\},Notoryctes_typhlops\{Test\})\{Test\})\{Test\})\{Test\},Monodelphis_domestica\{Test\},((Homo_sapiens,Mus_musculus),Loxodonta_africana));\
EOF\
\
# 3. Execute RELAX (numbat branch), RELAX (dietary guilds), and MEME (site-specific selection)\
for gene in TAS1R1 TAS1R2 TAS1R3; do\
    hyphy relax \\\
        --alignment "$HYPHY_DIR/$\{gene\}_18taxa_hyphy.fna" \\\
        --tree "$HYPHY_DIR/species_tree_numbat_test.nwk" \\\
        --test Test --reference Reference --models Minimal \\\
        --output "$HYPHY_DIR/$\{gene\}_numbat_relax.json"\
\
    hyphy relax \\\
        --alignment "$HYPHY_DIR/$\{gene\}_18taxa_hyphy.fna" \\\
        --tree "$HYPHY_DIR/species_tree_dietary_guilds.nwk" \\\
        --test Test --reference Reference --models Minimal \\\
        --output "$HYPHY_DIR/$\{gene\}_diet_relax.json"\
\
    hyphy meme \\\
        --alignment "$HYPHY_DIR/$\{gene\}_18taxa_hyphy.fna" \\\
        --tree "$HYPHY_DIR/species_tree_marsupials_only.nwk" \\\
        --branches Test --pvalue 0.05 \\\
        --output "$HYPHY_DIR/$\{gene\}_marsupials_meme.json"\
done\
}