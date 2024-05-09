## nanoMD

nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A sites, genes, and isoforms, and new mRNA.

## Table of Contents
<!-- TOC -->

- [nanoMD](#nanomd)
- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [nanoMD modules](#nanomd-modules)
    - [nanomd gene](#nanomd-gene)
    - [Usage](#usage)
    - [nanomd isoform](#nanomd-isoform)
    - [Usage](#usage)
    - [nanomd m6A sites](#nanomd-m6a-sites)
    - [Usage](#usage)
    - [nanomd new mRNA](#nanomd-new-mrna)
    - [Usage](#usage)
- [Scripts](#scripts)
- [Docker](#docker)
- [Conda Environment](#conda-environment)
- [Cite nanoMD](#cite-nanomd)

<!-- /TOC -->

## Overview

nanoMD  is designed to work with the output of the ONT direct RNA sequencing platform and can be used to identify new m6A sites, genes, and isoforms, as well as detect and quantify new mRNA. 

It is recommended to use the docker or conda environment to run the pipeline.

## Requirements

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

## nanoMD modules

nanomd.py is the main module of nanoMD, which includes the following sub-modules:
- gene
- isoform
- detectm6A
- detectnewmRNA

### nanomd gene

### Usage

`nanomd.py gene -i sample -f genome.cdna.fa -e 0.005 -o gene`

### nanomd isoform

### Usage

`nanomd.py isoform -i sample -f genome.cdna.fa -e 0.005 -o isoform`

### nanomd m6A sites

### Usage

`nanomd.py detectm6A -i sample -f genome.cdna.fa -e 0.005 -o m6A`

### nanomd new mRNA

### Usage

`nanomd.py detectnewmRNA -i sample -f genome.cdna.fa -e 0.005 -o newmRNA`

## Scripts

We provide a set of standalone scripts for 5EU detection and quantification.

### detect5EU.py

This script detects 5' untranslated regions (5EU) from the ONT direct RNA sequencing data.

### Usage

`python detect5EU.py -i sample.fastq -o 5EU.bed`

## Docker

If the user has docker installed, the following command can be used to run the pipeline in a docker container:

```
docker run -v /path/to/data:/data -it nanomd/nanomd:latest /bin/bash
```

## Conda Environment

If the user has conda installed, the following command can be used to create a conda environment for nanoMD:

1. Install conda
2. Create a new conda environment: `conda create -n nanomd python=3.6`
3. Activate the environment: `conda activate nanomd`
4. Install the required packages: `conda install -c bioconda minimap2 samtools bedtools flair tombo mines`
5. Install the required python packages: `pip install pandas numpy scipy sklearn matplotlib seaborn pysam`
6. Clone the nanoMD repository: `git clone https://github.com/epibiotek/nanomd.git`
7. Run the pipeline: `python nanomd/nanomd.py gene -i sample -f genome.cdna.fa -e 0.005 -o gene`

## Cite nanoMD

If you use nanoMD in your research, please cite the following paper: