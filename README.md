
## nanoMD

nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to
synchronously analyze the changes in m6A sites, genes, and isoforms, and newRNA

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
nanoMD  is designed to work with the output of the ONT direct RNA sequencing platform and can be used to identify new m6A sites, genes, and isoforms, as well as detect and quantify new RNA species.

## nanoMD modules
nanomd.py is the main module of nanoMD, which includes the following sub-modules:
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

## Docker

## Conda Environment

## Cite nanoMD
