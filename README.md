# scFusion

scFusion is a computational pipeline for detecting gene fusions at single-cell resolution. scFusion works on Linux/Mac OS.

## Software Prerequisite

The software below should be in your PATH.

- [STAR](https://github.com/alexdobin/STAR) >= 2.7.2d (tested on 2.7.2d and 2.7.8a)
- samtools (tested on version 1.10)
- bedtools
- python 3
- R >= 3.5 (tested on 3.5.1, 3.6.0, 4.0.2)

- python module: tensorflow (tested on version 2.3.0)
- python module: keras (tested on version 2.4.3)
- python module: [pyensembl](https://github.com/openvax/pyensembl)
- python module: numpy
- python module: scipy
- python module: pysam


## Recommend 

- 64 GB memory or more for each task

- 8 CPU cores or more for each task 

## Minimum 

- 64 GB Memory 

- 4 CPU cores

## Optional

- Job schedular (e.g. Slurm)

## Data Requirement

- Single cell RNA sequencing files, named with numbers (better with consecutive numbers). (*_1.fastq, *_2.fastq) (e.g. 1_1.fastq, 1_2.fastq, 2_1.fastq, 2_2.fastq) **The filenames must be \*_1.fastq and \*_2.fastq. \*_1.fq is not allowed.**

- STAR reference dataset(Please build it, see [Section 2 of STARManual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf))

- Reference genome file (*.fa)(like hg19.fa, file size = ~3G)

- GTF annotation file (*.gtf) (Can be obtained from Ensembl (ftp://ftp.ensembl.org/pub/), NCBI, or UCSC)

- [Mappabilityfile](https://genome.ucsc.edu/cgi-bin/hgTables) (Can be obtained from UCSC) **If you are using the hg38 version of mappability file, please ensure the file format is the same with hg19's. (Only 4 columns. chr, start, end, value)**

## Usage

First, unzip the hg19mappability file in the data folder.

Two versions of scFusion are included. The normal version (scFusion.py) runs the whole pipeline, while the job schedular version gives you a series of commands that you can run them with your job schedular's configuration. We recommend you using the job schedular version, since it can make use of all the available computational resources. 

Example:

If 300 pairs of files (501_1.fastq, 501_2.fastq, 502_1.fastq, 502_2.fastq, ..., 800_1.fastq, 800_2.fastq) in the testdata/, and you want to save the results at testout/, and the STAR reference folder is hg19StarIndex/, and 20 cores are available, please run:

`python software/scFusion.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex/ -t 20`

or the job schedular version:

`python software/scFusion_js.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex/ -t 20`

The results are shown on the folder: testout/FinalResult/FinalOutput*

### Expected running time

Running the test data of 10 cells costs about 10 minutes on an 8-core computer.

## Command Line Options

Parameters must be specified: 

    -f, --FileDir: The folder of single cell sequencing files
    
    -b, --Begin: The first index of single cell sequencing file
    
    -e, --End: The last index of single cell sequencing file
    
    -s, --STARReference: The reference folder of STAR. The reference should be built before running scFusion. 
    
Parameters with default values, but should be changed for your setting: 

    -g, --Genome: The genome reference file (*.fasta or *.fa), default is 'hg19.fa' in the data folder
    
    -a, --Annotation: The gtf annotation file (*.gtf), default is 'ref_annot.gtf'in the data folder
    
    -m, --Mappability: The mappability file, default is 'hg19mappability75.txt' in the data folder, left blank for keeping all reads
    
    -o, --OutDir: The output folder of the results and temperal files, default is the same as FileDir
    
    -t, --Thread: Number of threads can be used, at least 4, default is 8. It tells the total available threads that scFusion can apply.
    
    -l, --LimitThread: Number of maximum threads allowed for each STAR mapping task, default is 20
    
    -w, --Weight: The weight file of the deep-learning network, default is 'weight-V9-2.hdf5' in the data folder. If retraining is allowed, this file is the initial weight file of the network to be retrained; if retraining step is skipped, this weight file is used in the predicting step.  
    
    -E, --Epoch: The number of epochs in the retraining step between 3 and 999, default is 100
    
    -p, --Prefix: The prefix of result file, default is blank. This should be specified if users want to compare results of different settings.
    
    -v, --PvalueCutoff: Pvalue(FDR) cutoff of the statistical model, default is 0.05
    
    -n, --NetworkCutoff: Network score cutoff, default is 0.75
    
Step Controls:

    --SkipMapping: Skip STAR Mapping, if you already have the mapping result at OutDir/StarMapping/
    
    --SkipBS: Skip the basic processing step, if you already have *_FusionSupport.txt at OutDir/ChimericOut/
    
    --SkipCombining: Skip the combining step, if you already have *ChiDist_middle.txt at OutDir/ChiDist/
    
    --SkipRetrain: Skip the retraining step, and apply the weights specified by -w to the network
    
    --SkipPredict: Skip the predicting step using deep-learning network, if you already have the *ChiDist_filtered.txt at OutDir/ChiDist/


## If you are using job schedular

Take Slurm as example. The simplest execution command is `srun 1 20 XXXXX`.

First, run the scFusion_js.py

`python scFusion_js.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex_2.7.2b_normal/ -t 200`

Then, several commands are shown on the screen. Follow the instruction to run these commands. (one by one or parallelly)

if the command is `sh scFusion_dev/bin/CombinePipeline_before_FS.sh testout/ 601 605 scFusion_dev/bin/../data/ref_annot.gtf scFusion_dev/bin/../data/hg19mappability75.txt scFusion_dev/bin/../data/exon_probe.hg19.gene.new.bed scFusion_dev/bin/`, 

then run `srun 1 20 sh scFusion_dev/bin/CombinePipeline_before_FS.sh testout/ 601 605 scFusion_dev/bin/../data/ref_annot.gtf scFusion_dev/bin/../data/hg19mappability75.txt scFusion_dev/bin/../data/exon_probe.hg19.gene.new.bed scFusion_dev/bin/`


## Note
For non-academic use, please email Prof. Xi (ruibinxi@math.pku.edu.cn)
