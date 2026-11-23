{\rtf1\ansi\ansicpg1252\cocoartf2908
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Convergent pseudogenization of *TAS1R1* and divergent dietary selection across marsupial type 1 taste receptors\
\
This repository contains the custom Python and Bash scripts used to identify, splice, align, and perform phylogenetic and selection (*dN/dS*) analyses on *TAS1R1*, *TAS1R2*, and *TAS1R3* across 15 marsupial genomes and 3 placental outgroups.\
\
## Pipeline Overview\
1. **`01_extract_and_splice_tas1r.py`**: Runs `tblastn` against reference genome assemblies, clusters exon HSPs into genomic loci, extracts flanking genomic sequences, and performs splice-aware CDS reconstruction using `Exonerate`.\
2. **`02_align_and_build_trees.sh`**: Translates reconstructed CDS, aligns amino acid sequences with `MAFFT`, threads nucleotide codons back into frame-preserved alignments, and infers maximum-likelihood trees in `IQ-TREE`.\
3. **`03_verify_numbat_pseudogene.py`**: Audits all 45 marsupial coding sequences for frameshifts and premature stop codons, and extracts 61-bp genomic probes spanning the two nonsense mutations in *Myrmecobius fasciatus* *TAS1R1* for raw SRA read verification.\
4. **`04_run_hyphy_selection.sh`**: Masks stop codons, prepares species trees with branch labels, and executes `HyPhy RELAX` (numbat foreground branch and dietary-guild comparisons) and `HyPhy MEME` (site-specific episodic positive selection).\
5. **`05_plot_figures.py`**: Generates publication-ready vector (`PDF`, `SVG`) and high-resolution (`600 DPI PNG`) composite figures (**Fig. 1** and **Fig. 2**).\
}