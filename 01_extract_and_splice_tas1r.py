{\rtf1\ansi\ansicpg1252\cocoartf2908
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/env python3\
"""\
Step 1: Identify candidate TAS1R1, TAS1R2, and TAS1R3 genomic loci via tblastn\
and reconstruct full-length six-exon coding sequences using Exonerate.\
"""\
\
import os\
import glob\
import subprocess\
import tempfile\
from collections import defaultdict\
from Bio import SeqIO\
from Bio.Seq import Seq\
\
BASE_DIR = os.path.expanduser("~/tas1r_pipeline")\
GENOME_DIR = os.path.join(BASE_DIR, "genomes")\
QUERY_FAA = os.path.join(BASE_DIR, "queries/tas1r_monodelphis.faa")\
OUT_CDS_DIR = os.path.join(BASE_DIR, "results/cds")\
os.makedirs(OUT_CDS_DIR, exist_ok=True)\
\
LOCI_FNA = os.path.join(OUT_CDS_DIR, "tas1r_genomic_loci.fna")\
ALL_CDS_FNA = os.path.join(OUT_CDS_DIR, "tas1r_marsupials_cds.fna")\
REPORT_TSV = os.path.join(OUT_CDS_DIR, "tas1r_marsupials_final_report.tsv")\
\
# Load reference queries (Monodelphis domestica TAS1R1, TAS1R2, TAS1R3)\
queries = \{\}\
for rec in SeqIO.parse(QUERY_FAA, "fasta"):\
    for g in ("TAS1R1", "TAS1R2", "TAS1R3"):\
        if g in rec.description:\
            queries[g] = rec\
\
# 1. Run tblastn across all genome assemblies and cluster HSPs within 40 kb\
loci_records = []\
for fna_path in sorted(glob.glob(os.path.join(GENOME_DIR, "*.fna"))):\
    species = os.path.basename(fna_path).replace(".fna", "")\
    db_prefix = fna_path + "_db"\
    if not os.path.exists(db_prefix + ".nin"):\
        subprocess.run(["makeblastdb", "-in", fna_path, "-dbtype", "nucl", "-out", db_prefix], check=True)\
\
    cmd = [\
        "tblastn", "-query", QUERY_FAA, "-db", db_prefix,\
        "-evalue", "1e-10", "-outfmt", "6 qseqid sseqid pident length qlen sstart send bitscore"\
    ]\
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)\
\
    hits_by_gene = defaultdict(list)\
    for line in res.stdout.strip().splitlines():\
        qseqid, sseqid, pident, length, qlen, sstart, send, bitscore = line.split("\\t")\
        gene = next((g for g in ("TAS1R1", "TAS1R2", "TAS1R3") if g in qseqid), None)\
        if gene and float(pident) >= 55.0:\
            s0, s1 = int(sstart), int(send)\
            strand = "+" if s0 < s1 else "-"\
            hits_by_gene[gene].append((sseqid, strand, min(s0, s1), max(s0, s1), float(bitscore)))\
\
    genome_dict = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))\
    for gene, hsps in hits_by_gene.items():\
        # Group HSPs by scaffold and strand, pick highest-scoring cluster within 40 kb\
        clusters = defaultdict(list)\
        for scaf, strand, start, end, score in hsps:\
            clusters[(scaf, strand)].append((start, end, score))\
\
        best_key = max(clusters, key=lambda k: sum(x[2] for x in clusters[k]))\
        scaf, strand = best_key\
        coords = clusters[best_key]\
        loc_start = max(0, min(x[0] for x in coords) - 1500)\
        loc_end = min(len(genome_dict[scaf].seq), max(x[1] for x in coords) + 1500)\
\
        subseq = genome_dict[scaf].seq[loc_start:loc_end]\
        if strand == "-":\
            subseq = subseq.reverse_complement()\
\
        header = f"\{species\}|\{gene\}|\{scaf\}:\{loc_start+1\}-\{loc_end\}(\{strand\})"\
        loci_records.append(f">\{header\}\\n\{str(subseq)\}\\n")\
\
with open(LOCI_FNA, "w") as out:\
    out.writelines(loci_records)\
\
# 2. Reconstruct spliced CDS via Exonerate protein2genome\
cds_records = []\
report_rows = ["Species\\tGene\\tScaffold_Range\\tStatus\\tCDS_bp\\tAA_len\\tInternal_Stops\\n"]\
\
for loc_rec in SeqIO.parse(LOCI_FNA, "fasta"):\
    species, gene, scaf_range = loc_rec.id.split("|")\
    ref_rec = queries[gene]\
\
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as qf, \\\
         tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:\
        SeqIO.write(ref_rec, qf, "fasta")\
        SeqIO.write(loc_rec, tf, "fasta")\
        q_path, t_path = qf.name, tf.name\
\
    exo = subprocess.run([\
        "exonerate", "--model", "protein2genome", "--bestn", "1",\
        "--showalignment", "no", "--showvulgar", "no",\
        "--ryo", ">%ti\\n%tcs\\n", q_path, t_path\
    ], capture_output=True, text=True, check=True)\
    os.remove(q_path)\
    os.remove(t_path)\
\
    nt_lines = [l.strip() for l in exo.stdout.splitlines() if l and not l.startswith(">") and not l.startswith("Command") and not l.startswith("Hostname") and not l.startswith("--")]\
    spliced_nt = "".join(nt_lines).upper()\
    trim_len = len(spliced_nt) - (len(spliced_nt) % 3)\
    aa_seq = str(Seq(spliced_nt[:trim_len]).translate(to_stop=False))\
    body_aa = aa_seq[:-1] if aa_seq.endswith("*") else aa_seq\
    stops = [i + 1 for i, c in enumerate(body_aa) if c == "*"]\
\
    status = "PSEUDOGENE" if stops else "INTACT"\
    cds_records.append(f">\{species\}|\{gene\}|status=\{status\}|aa_len=\{len(body_aa)\}|stops=\{len(stops)\}\\n\{spliced_nt\}\\n")\
    report_rows.append(f"\{species\}\\t\{gene\}\\t\{scaf_range\}\\t\{status\}\\t\{len(spliced_nt)\}\\t\{len(body_aa)\}\\t\{stops\}\\n")\
\
with open(ALL_CDS_FNA, "w") as out:\
    out.writelines(cds_records)\
with open(REPORT_TSV, "w") as out:\
    out.writelines(report_rows)\
\
print("Completed CDS reconstruction and audit report.")\
}