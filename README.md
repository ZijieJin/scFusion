# scFusion  [![DOI](https://zenodo.org/badge/372129480.svg)](https://zenodo.org/badge/latestdoi/372129480)

scFusion is a computational pipeline for detecting gene fusions at single-cell resolution. scFusion works on Linux/Mac OS. If you have any questions related to scFusion, please visit https://github.com/ZijieJin/scFusion and post them on the *Issues* page.

## Software Prerequisite

The software below should be in your PATH. **(They can all be installed by conda and pip)**

- [STAR](https://github.com/alexdobin/STAR) >= 2.7.2d (tested on 2.7.2d and 2.7.8a)
- samtools (tested on version 1.10)
- bedtools (tested on version 2.29.2)
- python 3
- R >= 3.5 (tested on 3.5.1, 3.6.0, 4.0.2)
- R package: stringr
- python module: tensorflow (tested on version 2.3.0)
- python module: keras (tested on version 2.4.3)
- python module: [pyensembl](https://github.com/openvax/pyensembl) (tested on version 1.8.8)
- python module: numpy (tested on version 1.18.5)
- python module: scipy (tested on version 1.4.1 and 1.5.0)
- python module: pysam (tested on version 0.16.0.1)


## Recommend 

- 64 GB memory or more for each task

- 8 CPU cores or more for each task 

## Minimum 

- 64 GB Memory 

- 4 CPU cores

## Optional

- Job scheduler (e.g. Slurm)

## Data Requirement

- Single cell RNA sequencing files, named with numbers (better with consecutive numbers). (*_1.fastq, *_2.fastq) (e.g. 1_1.fastq, 1_2.fastq, 2_1.fastq, 2_2.fastq) **The filenames must be \*_1.fastq and \*_2.fastq. \*_1.fq is not allowed.**

- STAR reference dataset(Please build it, see [Section 2 of STARManual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf))

- Reference genome file (*.fa)(like hg19.fa, file size = ~3G)

- GTF annotation file (*.gtf) (Can be obtained from Ensembl (ftp://ftp.ensembl.org/pub/), NCBI, or UCSC)

- [Mappabilityfile](https://genome.ucsc.edu/cgi-bin/hgTables) (Can be obtained from UCSC) **If you are using the hg38 version of mappability file, please ensure the file format is the same as hg19's. (Only 4 columns. chr, start, end, value)**

## Usage

First, unzip the hg19mappability file in the data folder.

Two versions of scFusion are included. The normal version (scFusion.py) runs the whole pipeline, while the job schedular version gives you a series of commands that you can run with your job schedular's configuration. We recommend you using the job scheduler version since it can make use of all the available computational resources. 

Example:

If 300 pairs of files (501_1.fastq, 501_2.fastq, 502_1.fastq, 502_2.fastq, ..., 800_1.fastq, 800_2.fastq) in the testdata/, and you want to save the results at testout/, and the STAR reference folder is hg19STARIndex/, and 20 cores are available, please run:

`python software/scFusion.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19STARIndex/ -t 20 -g hg19.fa -a ref_annot.gtf`

or the job schedular version:

`python software/scFusion_js.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19STARIndex/ -t 20 -g hg19.fa -a ref_annot.gtf`

The results are shown on the folder: testout/FinalOutput.abridged.txt and testout/FinalOutput.full.txt

## The test dataset

The test dataset helps you to check whether you have correctly installed scFusion. The testdata is only for installation checking, it does not have any biological meanings.

### Command to run

Please modify the command below with parameters fitting your environment.

`python software/scFusion.py -f testdata/ -o testout/ -b 1 -e 10 -s hg19STARIndex/ -t 8 -n 0.9 -g hg19.fa -a ref_annot.gtf`

### Expected running time and result

Running the test data of 10 cells costs about 30 minutes (including the one-time file preparing time) on an 8-core computer.

scFusion is expected to report IGHJ5-IGHA1 fusion in the testout/FinalResult/Final.txt file.

## Command Line Options

Parameters must be specified: 

    -f, --FileDir: The folder of single-cell sequencing files
    
    -b, --Begin: The first index of single-cell sequencing file
    
    -e, --End: The last index of single-cell sequencing file
    
    -s, --STARReference: The reference folder of STAR. The reference should be built before running scFusion. 
    
Parameters with default values, but should be changed for your setting: 

    -g, --Genome: The genome reference file (*.fasta or *.fa), default is 'hg19.fa' in the data folder
    
    -a, --Annotation: The gtf annotation file (*.gtf), default is 'ref_annot.gtf'in the data folder
    
    -m, --Mappability: The mappability file, default is 'hg19mappability75.txt' in the data folder, left blank for keeping all reads
    
    -o, --OutDir: The output folder of the results and temporal files, default is the same as FileDir
    
    -t, --Thread: Number of threads can be used, at least 4, default is 8. It tells the total available threads that scFusion can apply.
    
    -l, --LimitThread: Number of maximum threads allowed for each STAR mapping task, default is 20
    
    -w, --Weight: The weight file of the deep-learning network, default is 'weight-V9-2.hdf5' in the data folder. If retraining is allowed, this file is the initial weight file of the network to be retrained; if the retraining step is skipped, this weight file is used in the predicting step.  
    
    -E, --Epoch: The number of epochs in the retraining step between 3 and 999, default is 10
    
    -p, --Prefix: The prefix of result file, default is blank. This should be specified if users want to compare the results of different settings.
    
    -v, --PvalueCutoff: Pvalue(FDR) cutoff of the statistical model, default is 0.05
    
    -n, --NetworkCutoff: Network score cutoff, default is 0.75
    
Step Controls:
    
    --Rename: Rename the files with consecutive indexes
    
    --SkipMapping: Skip STAR Mapping, if you already have the mapping result at OutDir/StarMapping/
    
    --SkipBS: Skip the basic processing step, if you already have *_FusionSupport.txt at OutDir/ChimericOut/
    
    --SkipCombining: Skip the combining step, if you already have *ChiDist_middle.txt at OutDir/ChiDist/
    
    --SkipRetrain: Skip the retraining step, and apply the weights specified by -w to the network
    
    --SkipPredict: Skip the predicting step using the deep-learning network, if you already have the *ChiDist_filtered.txt at OutDir/ChiDist/


## If you are using a job scheduler

Take Slurm as an example. The simplest execution command is `srun 1 20 XXXXX`.

First, run the scFusion_js.py

`python scFusion_js.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex_2.7.2b_normal/ -t 200`

Then, several commands are shown on the screen. Follow the instruction to run these commands. (one by one or parallelly)

if the command is `sh scFusion_dev/bin/CombinePipeline_before_FS.sh testout/ 601 605 scFusion_dev/bin/../data/ref_annot.gtf scFusion_dev/bin/../data/hg19mappability75.txt scFusion_dev/bin/../data/exon_probe.hg19.gene.new.bed scFusion_dev/bin/`, 

then run `srun 1 20 sh scFusion_dev/bin/CombinePipeline_before_FS.sh testout/ 601 605 scFusion_dev/bin/../data/ref_annot.gtf scFusion_dev/bin/../data/hg19mappability75.txt scFusion_dev/bin/../data/exon_probe.hg19.gene.new.bed scFusion_dev/bin/`


## Alternative installation

The above installation is easy and we use Python, R, and Shell to make scFusion easy to be used. All prerequisites can be installed by conda and pip. If you have trouble with the installation, we provide a Docker image that contains all software pre-installed for running scFusion. it is available here: https://hub.docker.com/r/jzj2035198/scfusion . 

If you have docker installed, you can pull the image like so:

`docker pull jzj2035198/scfusion`

Or you can build the image using Docker/Dockerfile.

Below we assume all the required files and folders are in the directory XXX, run

`docker run -v XXX:/data --rm jzj2035198/scfusion python -u /usr/local/src/scFusion/scFusion.py -f /data/testdata/ -o /data/testout/ -b 501 -e 800 -s /data/hg19StarIndex/ -t 20 -g /data/hg19.fa -a /data/ref_annot.gtf`

The `XXX:/data` means you map the XXX folder to /data, so all your files and directories in XXX can be found in /data. 

### To test docker

run `docker run -v XXX:/data --rm jzj2035198/scfusion python -u /usr/local/src/scFusion/scFusion.py -f /usr/local/src/scFusion/Testdata/ -o /data/out/ -b 1 -e 10 -s /data/hg19StarIndex/ -t 20 -n 0.9`

The result file locates at XXX/FinalOutput.abridged.txt or XXX/FinalResult/FinalOutput.abridged.txt

## About the annotation file

The annotation file (\*.gtf) may have different format, so making scFusion be compatible with all the formats is difficult. 'gene_name' is the gene name indicator in the annotation, and 'gene_type' or 'gene_biotype' are the gene type (pseudo gene or LncRNA). The annotation file containing only chr1-chrY can work best with scFusion.

## Commercial usage

For non-academic use, please email Prof. Xi (ruibinxi@math.pku.edu.cn) to obtain the paid commercial license.

For academic use, source code is licensed under MIT License. 
