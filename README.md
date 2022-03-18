# scFusion  [![DOI](https://zenodo.org/badge/372129480.svg)](https://zenodo.org/badge/latestdoi/372129480)

scFusion is a computational pipeline for detecting gene fusions at single-cell resolution. scFusion works on Linux/Mac OS. If you have any questions related to scFusion, please visit https://github.com/ZijieJin/scFusion and post them on the *Issues* page or email me: jinzijie@pku.edu.cn

## Software Prerequisite

The software below should be in your PATH. **(They can be installed by conda and pip)**

- [STAR](https://github.com/alexdobin/STAR) >= 2.7.2d (tested on 2.7.2d and 2.7.8a)
- samtools (tested on version 1.10)
- bedtools (tested on version 2.29.2)
- python 3
- R >= 3.5 (tested on 3.5.1, 3.6.0, 4.0.2, 4.0.5)
- R package: stringr
- python module: tensorflow (tested on version 2.3.0)
- python module: keras (tested on version 2.4.3)
- python module: [pyensembl](https://github.com/openvax/pyensembl) (tested on version 1.8.8)
- python module: numpy (tested on version 1.18.5)
- python module: scipy (tested on version 1.4.1 and 1.5.0)
- python module: pysam (tested on version 0.16.0.1)


## Recommend Configuration

- 64 GB memory or more for each task

- 8 CPU cores or more for each task 

## Minimum Configuration

- 64 GB Memory 

- 4 CPU cores

## Optional Configuration

- Job scheduler (e.g. Slurm)

## Data Requirement

- Single cell RNA sequencing files from Smart-Seq protocol. File names should be *_1.fastq, *_2.fastq (e.g. 1_1.fastq, 1_2.fastq, 2_1.fastq, 2_2.fastq) ** \*_1.fq is not allowed.**

- Reference genome file (*.fa)(like hg19.fa, file size = ~3G)

- GTF annotation file (*.gtf) (Can be obtained from Ensembl (ftp://ftp.ensembl.org/pub/), NCBI, or UCSC)

- (Optional) [Mappabilityfile](https://genome.ucsc.edu/cgi-bin/hgTables) (Can be obtained from UCSC) **If you are using the hg38 version of mappability file, please ensure the file format is the same as hg19's. (Only 4 columns. chr, start, end, value)** If you do not provide this file, scFusion will turn off the mappability filter.

## Usage

Download all scripts, unzip the hg19mappability file in the data folder, and run scFusion.py in the command line.

See **Manual.pdf** for details (it will be uploaded soon).


## Alternative installation

scFusion is easy to use, consisting of Python, R, and Shell scripts. All prerequisites can be installed by conda and pip. If you have trouble with the installation, we provide a Dockerfile to build a Docker image. 

Below we assume all the required files and folders are in the directory XXX, run

`docker run -v XXX:/data --rm jzj2035198/scfusion python -u /usr/local/src/scFusion/scFusion.py [commands]`

The `XXX:/data` means you map the XXX folder to /data, so all your files and directories in XXX can be found in /data. 

## About the annotation file

The annotation file (\*.gtf) may have different format, so making scFusion be compatible with all the formats is difficult. 'gene_name' is the gene name indicator in the annotation, and 'gene_type' or 'gene_biotype' are the gene type (pseudo gene or LncRNA). 

## Commercial use

For non-academic use, please email Prof. Xi (ruibinxi@math.pku.edu.cn) to obtain the paid commercial license.

For academic use, source code is licensed under MIT License. 
