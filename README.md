# scFusion

scFusion is a computational pipeline for detecting gene fusions at single-cell resolution. 

## Software Prerequisite

The software below should be in your PATH.

- [STAR](https://github.com/alexdobin/STAR) >= 2.7 
- samtools (tested with version 0.1.19, and 1.2 should work)
- bedtools
- python 3
- R >= 3.5

- python module: tensorflow
- python module: keras
- python module: [pyensembl](https://github.com/openvax/pyensembl)
- python module: numpy
- python module: scipy
- python module: pysam


## Recommand

- Job scheduler (e.g. Slurm)

- 64 GB Memory or more for each task

- 200 CPU cores or more

## Minimum

- 64 GB Memory

- 20 CPU cores

## Data Requirement

- A Series of single cell sequencing file, renamed with numbers, recommand with continuous numbers. (*_1.fastq, *_2.fastq) (like 1_1.fastq, 1_2.fastq, 2_1.fastq, 2_2.fastq)

- STAR reference dataset(Please build it, see [Section 2 of STARManual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf))

- Reference genome file (*.fa)(like hg19.fa, file size = ~3G)

- GTF annotation file (*.gtf) (Can be obtained from Ensembl (ftp://ftp.ensembl.org/pub/), NCBI, or UCSC)

- [Mappabilityfile](https://genome.ucsc.edu/cgi-bin/hgTables) (Can be obtained from UCSC)

## Usage

Uncompress the badgene.zip to the data/badgene/ folder. 

Two versions of scFusion are contained. The normal version (scFusion.py) run the whole pipelines, while the job scheduler version give you a series of commands that you can run them with your job scheduler's configuration. We recommand you to use the job scheduler version, since it can make use of all the computational resources available. 

Example:

if 300 sequencing data file (501_1.fastq, 501_2.fastq, 502_1.fastq, 502_2.fastq, ..., 800_1.fastq, 800_2.fastq) in the testdata/, and you want to save the results at testout/, and the STAR reference folder is hg19StarIndex_2.7.2b_normal/, and 200 cores are available, please run:

`python scFusion.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex_2.7.2b_normal/ -t 200`

or

`python scFusion_js.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex_2.7.2b_normal/ -t 200`

## Command Line Options

Parameters must be specified: 

    -f, --FileDir: The folder of single cell sequencing files
    
    -b, --Begin: The first index of single cell sequencing file
    
    -e, --End: The last index of single cell sequencing file
    
    -s, --STARReference: The reference folder of STAR. The reference should be built before running scFusion. 
    
Parameters with default value, but should be change for your setting: 

    -g, --Genome: The genome reference file (*.fasta or *.fa), default is the 'hg19.fa' in the data folder
    
    -a, --Annotation: The gtf annotation file (*.gtf), default is the 'ref_annot.gtf'in the data folder
    
    -m, --Mappability: The mappability file, default is the 'hg19mappability75.txt' in the data folder, left blank for keep all reads
    
    -o, --OutDir: The output folder of the results and temperal files, default is the same as FileDir
    
    -t, --Thread: Number of threads can be used, at least 4, default is 8. It tells the total available threads that scFusion can apply.
    
    -l, --LimitThread: Number of maximum threads allowed for each STAR mapping task, default is 20
    
    -w, --Weight: The weight file of network, default is the 'weight-V9-2.hdf5' in the data folder. If Retraining is allowed, the file is the initial weight of the network to be retrained; if Retraining step is skipped, this weight file is used in the predict.  
    
    -E, --Epoch: The number of epoch in the retraining step between 3 and 999, default is 100
    
    -p, --Prefix: The prefix of result file, default is blank. This should be specified if users want to compare results of diffrerent settings.
    
    -v, --PvalueCutoff: Pvalue(FDR) cutoff of the statistical model, default is 0.05
    
    -N, --NetworkCutoff: Network classification probability cutoff, default is 0.75
    
Step Controls:

    --SkipMapping: Skip STAR Mapping, if you already have the mapping result at OutDir/StarMapping/
    
    --SkipBS: Skip the basic processing step, if you already have the *_FusionSupport.txt at OutDir/ChimericOut/
    
    --SkipCombining: Skip the combining step, if you already have the *ChiDist_middle.txt at OutDir/ChiDist/
    
    --SkipRetrain: Skip the retraining step, and apply the weights specified by -w to the network
    
    --SkipPredict: Skip the predicting step using network, if you already have the *ChiDist_filtered.txt at OutDir/ChiDist/


## If you are using job scheduler

Take Slurm as example. The simplest execution command is `srun 1 20 XXXXX`.

First, run the scFusion_js.py

`python scFusion_js.py -f testdata/ -o testout/ -b 501 -e 800 -s hg19StarIndex_2.7.2b_normal/ -t 200`

Then, several commands are shown on the screen. Follow the instruction to run these commands. (one by one or parallelly)

if the command is `sh scFusion_dev/bin/CombinePipeline_before_FS.sh testout/ 601 605 scFusion_dev/bin/../data/ref_annot.gtf scFusion_dev/bin/../data/hg19mappability75.txt scFusion_dev/bin/../data/exon_probe.hg19.gene.new.bed scFusion_dev/bin/`, 

then run `srun 1 20 sh scFusion_dev/bin/CombinePipeline_before_FS.sh testout/ 601 605 scFusion_dev/bin/../data/ref_annot.gtf scFusion_dev/bin/../data/hg19mappability75.txt scFusion_dev/bin/../data/exon_probe.hg19.gene.new.bed scFusion_dev/bin/`


## Note
For non-academic use, please email Prof. Xi (ruibinxi@math.pku.edu.cn)
