# nanoMD

nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A sites, genes, and isoforms, and new mRNA.

# Table of Contents
<!-- TOC -->

- [nanoMD](#nanomd)
- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [Requirements](#requirements)
- [nanoMD modules](#nanomd-modules)
    - [nanomd gene](#nanomd-gene)
    - [Usage](#usage)
    - [nanomd isoform](#nanomd-isoform)
    - [Usage](#usage)
    - [nanomd modification sites](#nanomd-modification-sites)
    - [Usage](#usage)
    - [nanomd nascentRNA](#nanomd-nascentrna)
    - [Usage](#usage)
- [Scripts](#scripts)
    - [detect5EU.py](#detect5eupy)
    - [Usage](#usage)
- [Docker](#docker)
- [Conda Environment](#conda-environment)
- [Cite nanoMD](#cite-nanomd)

<!-- /TOC -->

# Overview

nanoMD  is designed to work with the output of the ONT direct RNA sequencing platform and can be used to identify new m6A sites, genes, and isoforms, as well as detect and quantify new mRNA. 

It is recommended to use the docker or conda environment to run the pipeline.

# Requirements

1. Python 3.6+
2. Python modules:
    - pandas
    - numpy
    - scipy
    - sklearn
    - matplotlib
    - seaborn
    - pysam
3. minimap2
4. samtools
5. bedtools
6. flair
7. tombo
8. mines
9. metaPlotR

# nanoMD modules

nanomd.py is the main module of nanoMD, which includes the following sub-modules:
- gene
- isoform
- detectMod
- nascentRNA

## nanomd gene

## Usage

`nanomd gene -i ../input/{}/pass.fq.gz -r ../reference/fasta/genome.fa -o ./{}_gene.sam --parms '--secondary=no --cs -a --sam-hit-only'`

## nanomd isoform

## Usage

`nanomd isoform -i ../input/{}/pass.fq.gz -r ../reference/fasta/transcripts.fa -o ./{}_transcripts.sam`

## nanomd modification sites

## Usage

`nanomd detectMod -i ../input/{}/pass.fq.gz -s ../01_map_gene/{}_gene.sam -b ../reference/genes/genes.bed -r ../reference/genes/region_sizes.txt  -p {}`

## nanomd nascentRNA

## Usage

`nanomd nascentRNA -i ../input/{}/pass.fq.gz -s ../01_map_gene/{}_gene.sam -b U -m ~/soft/newRNA.pkl -p {}`

# Scripts

We provide a set of standalone scripts for 5EU detection and quantification.

## detect5EU.py

This script detects 5' untranslated regions (5EU) from the ONT direct RNA sequencing data.

## Usage

`python detect5EU.py -i sample.fastq -o 5EU.bed`

# Docker

If the user has docker installed, the following command can be used to run the pipeline in a docker container:

```
docker run -v /path/to/data:/data -it nanomd/nanomd:latest /bin/bash
```

# Conda Environment

If the user has conda installed, the following command can be used to create a conda environment for nanoMD:

1. Install conda
2. Create a new conda environment: `conda create -n nanomd python=3.6`
3. Activate the environment: `conda activate nanomd`
4. Install the required packages: `conda install -c bioconda minimap2 samtools bedtools flair tombo mines`
5. Install the required python packages: `pip install pandas numpy scipy sklearn matplotlib seaborn pysam`
6. Clone the nanoMD repository: `git clone https://github.com/epibiotek/nanomd.git`
7. Run the pipeline: `python nanomd/nanomd.py gene -i sample -f genome.cdna.fa -e 0.005 -o gene`

# Cite nanoMD

If you use nanoMD in your research, please cite the following paper: