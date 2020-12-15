from __future__ import print_function
from __future__ import division
import sys
import getopt
import os
import threading
import numpy
import subprocess


def help():
    print("\nThe job scheduler version of scFusion gives you the running commands,\n and you should apply them with your own job scheduler's configurations\n")
    print('----------Command Line Options: ----------')
    print('Parameters must be specified: ')
    print('-f, --FileDir: The folder of single cell sequencing files')
    print('-b, --Begin: The first index of single cell sequencing file')
    print('-e, --End: The last index of single cell sequencing file')
    print('-s, --STARReference: The reference folder of STAR')
    print('Parameters with default value, but should be change for your setting: ')
    print('-g, --Genome: The genome reference file (*.fasta or *.fa), default is in the data folder')
    print('-a, --Annotation: The gtf annotation file (*.gtf), default is in the data folder')
    print('-m, --Mappability: The mappability file, default is in the data folder, left blank for keep all reads')
    print('-o, --OutDir: The output folder of the results, default is the same as FileDir')
    print('-t, --Thread: Number of thread, at least 4, default is 8')
    print('-l, --LimitThread: Number of maximum threads allowed for each task, default is 20')
    print('-w, --Weight: The weight file of network, default is in the data folder')
    print('-E, --Epoch: The number of epoch in the retraining step between 3 and 999, default is 100')
    print('-p, --Prefix: The prefix of result file, default is blank')
    print('-v, --PvalueCutoff: Pvalue cutoff, default is 0.05')
    print('-N, --NetworkCutoff: Network classification probability cutoff, default is 0.75')
    print('Step Controls:')
    print('--Rename: Rename the files with the continuous index')
    print('--SkipMapping: Skip STAR Mapping, if you already have the mapping result at OutDir/StarMapping/')
    print('--SkipBS: Skip the basic processing step')
    print('--SkipCombining: Skip the combining step')
    print('--SkipRetrain: Skip the retraining step, and apply the weights specified by -w to the network')
    print('--SkipPredict: Skip the predicting step using network')


def RunSTARMapping(s, e, core):
    print('sh ' + codedir + 'StarMapping_Chimeric.sh ' + filedir + ' ' + str(s) + ' ' + str(e) + ' ' + outdir +
              'StarMapping/ ' + starrefdir + ' ' + str(core) + '\n')


def RunBS(s, e):
    print('sh ' + codedir + 'CombinePipeline_before_FS.sh ' + outdir + ' ' + str(s) + ' ' + str(e) + ' ' +
              gtffilepath + ' ' + mappabilityfilepath + ' ' + exonposfilepath + ' ' + codedir + '\n')



try:
    lastd = sys.argv[0].rfind('/')
    if lastd == -1:
        codedir = './bin/'
    else:
        codedir = sys.argv[0][:lastd] + '/bin/'
    numthread = 8
    maxthread = 20
    filedir = ''
    outdir = ''
    starrefdir = ''
    SkipMapping = False
    SkipBS = False
    SkipCombining = False
    SkipRetrain = False
    SkipPredict = False
    Rename = False
    hg19filepath = codedir + '../data/hg19.fa'
    gtffilepath = codedir + '../data/ref_annot.gtf'
    mappabilityfilepath = codedir + '../data/hg19mappability75.txt'
    exonposfilepath = codedir + '../data/exon_probe.hg19.gene.new.bed'
    weightfilepath = codedir + '../data/weight-V9-2.hdf5'
    pv = 0.05
    FakeProb = 0.75
    epoch = 100
    prefix = '.'
    start = -10000000
    end = -10000000
    cellindex = []
    opts, args = getopt.getopt(sys.argv[1:], 'ht:s:g:a:b:e:f:m:p:w:v:n:o:l:E:',
                               ['help', 'Thread=', 'STARReference=', 'Genome=', 'Annotation=', 'Begin=', 'End=',
                                'FileDir=', 'Mappability=', 'LimitThread=',
                                'Prefix=', 'Weight=', 'Epoch=', 'PvalueCutoff=', 'NetworkCutoff=', 'OutDir=',
                                'SkipMapping',
                                'SkipBS', 'SkipCombining', 'SkipRetrain', 'SkipPredict', 'Rename'])
    for opt, arg in opts:
        if opt == '-h' or opt == '--help':
            help()
            sys.exit()
        elif opt in ('-t', '--Thread'):
            numthread = int(arg)
        elif opt in ('-l', '--LimitThread'):
            maxthread = int(arg)
        elif opt in ('-f', '--FileDir'):
            filedir = arg
        elif opt in ('-o', '--OutDir'):
            outdir = arg
        elif opt in ('-s', '--STARReference'):
            starrefdir = arg
        elif opt in ('-g', '--Genome'):
            hg19filepath = arg
        elif opt in ('-a', '--Annotation'):
            gtffilepath = arg
        elif opt in ('-m', '--Mappability'):
            mappabilityfilepath = arg
        elif opt in ('-w', '--Weight'):
            weightfilepath = arg
        elif opt in ('-b', '--Begin'):
            start = int(arg)
        elif opt in ('-e', '--End'):
            end = int(arg)
        elif opt in ('-p', '--Prefix'):
            prefix = arg
        elif opt in ('-v', '--PvalueCutoff'):
            pv = float(arg)
        elif opt in ('-n', '--NetworkCutoff'):
            FakeProb = float(arg)
        elif opt in ('-E', '--Epoch'):
            epoch = int(arg)
        elif opt == '--SkipMapping':
            SkipMapping = True
        elif opt == '--SkipBS':
            SkipBS = True
        elif opt == '--SkipCombining':
            SkipCombining = True
        elif opt == '--SkipRetrain':
            SkipRetrain = True
        elif opt == '--SkipPredict':
            SkipPredict = True
        elif opt == '--Rename':
            Rename = True
    # printpar()
    if numthread < 4:
        sys.stderr.write('Please set a valid number of threads\n')
        help()
        sys.exit()
    if maxthread < 1:
        sys.stderr.write('Please set a valid number of threads\n')
        help()
        sys.exit()
    if filedir == '' or not os.path.exists(filedir):
        sys.stderr.write('Please set the valid file folder path\n')
        help()
        sys.exit()
    if not 2 < epoch < 1000:
        sys.stderr.write('Please set a valid number of epoch\n')
        help()
        sys.exit()
    filedir += '/'
    if outdir == '':
        outdir = filedir
    os.system('mkdir -p ' + outdir)
    if starrefdir == '' or not os.path.exists(starrefdir):
        sys.stderr.write('Please set the valid STAR reference path\n')
        help()
        sys.exit()
    if hg19filepath == '' or not os.path.exists(hg19filepath):
        sys.stderr.write('Please set the valid Genome reference path\n')
        help()
        sys.exit()
    if gtffilepath == '' or not os.path.exists(gtffilepath):
        sys.stderr.write('Please set the valid annotation reference path\n')
        help()
        sys.exit()
    if mappabilityfilepath != '' and not os.path.exists(mappabilityfilepath):
        sys.stderr.write('Please set the valid mappability file path\n')
        help()
        sys.exit()
    if weightfilepath != '' and not os.path.exists(weightfilepath):
        sys.stderr.write('Please set the valid weight file path\n')
        help()
        sys.exit()
    if start < 0:
        sys.stderr.write('Please set the valid first index of single cells\n')
        help()
        sys.exit()
    if end < 0 or end < start:
        sys.stderr.write('Please set the valid last index of single cells\n')
        help()
        sys.exit()

    print(
        "\nThe job scheduler version of scFusion gives you the running commands,\n and you should apply them with your own job scheduler's configurations\n")
    print('\nPreparing for scFusion!\n')
    os.system('python ' + codedir + 'GetGenePos.py ' + gtffilepath + ' > ' + codedir + '/../data/GenePos.txt\n')
    os.system('python ' + codedir + 'GetExonPos.py ' + gtffilepath + ' > ' + exonposfilepath + '\n')
    aaa = subprocess.check_output(
        'pyensembl install --reference-name GRCH37 --annotation-name my_genome_features --gtf ' + gtffilepath, shell=True, stderr=subprocess.STDOUT)
    for i in range(start, end + 1):
        if os.path.exists(filedir + str(i) + '_2.fastq'):
            cellindex.append(i)
    numcell = len(cellindex)
    os.system('mkdir -p ' + outdir + '/StarMapping/\n')
    os.system('mkdir -p ' + outdir + '/ChimericOut/\n')
    os.system('mkdir -p ' + outdir + '/Expr/\n')
    print('\nParameter Check Complete!\n')
    if Rename:
        print('Please run this commmand to rename files!\n--------------------------------------------')
        print('python ' + codedir + 'RenameFastqFiles.py ')
        print('--------------------------------------------\n\n')
    if not SkipMapping:
        if numthread <= min(24, maxthread):
            print(
                'Please run this commmand with your job scheduler!\n--------------------------------------------')
            RunSTARMapping(start, end, numthread - 1)
            print('--------------------------------------------\n\n')
        else:
            numtask = min(round(numthread / min(16, maxthread)), numcell)
            numcore = int(min(numpy.floor((numthread - 1) / numtask), maxthread))
            numcelleachtask = int(numpy.ceil(numcell / numtask))
            print('Please run these commmands parallelly with your job scheduler! \n(Hint: a for loop is helpful '
                  'since these commands are the same except the cell index)\n--------------------------------------------')
            for j in range(numtask):
                RunSTARMapping(cellindex[j * numcelleachtask], min(cellindex[min((j + 1) * numcelleachtask - 1, numcell - 1)], end), numcore)
            print('--------------------------------------------\n\n')

    if not SkipBS:
        numtask = min(numthread - 1, numcell)
        numcelleachtask = int(numpy.ceil(numcell / numtask))
        print('Please run these commmands parallelly with your job scheduler! \n(Hint: a for loop is helpful since '
              'these commands are the same except the cell index)\n--------------------------------------------')
        for j in range(numtask):
            RunBS(cellindex[j * numcelleachtask], min(cellindex[(j + 1) * numcelleachtask - 1], end))
        print('--------------------------------------------\n\n')

    print('Please run the commands below one by one with your job scheduler!')
    print('--------------------------------------------')
    if not SkipCombining:
        print('sh ' + codedir + 'CombinePipeline_startwith_FS.sh ' + outdir + ' ' + str(start) + ' ' + str(end) + ' ' + prefix + ' ' +
            weightfilepath + ' ' + hg19filepath + ' ' + gtffilepath + ' ' + codedir + '\n')
    if not SkipRetrain:
        print('sh ' + codedir + 'CombinePipeline_Retrain.sh ' + outdir + ' ' + prefix + ' ' + weightfilepath + ' ' + str(epoch) + ' ' + codedir + '\n')
        print('Please check the output weights in the Retrain folder and specify the weight file you use in the next step!\n')
    if not SkipPredict:
        print('sh ' + codedir + 'CombinePipeline_Predict.sh ' + outdir + ' ' + str(start) + ' ' + str(
            end) + ' ' + prefix + ' YOUR_WEIGHT_FILE ' + hg19filepath + ' ' + gtffilepath + ' ' + codedir + '\n')
    print(
        'sh ' + codedir + 'CombinePipeline_startwith_ChiDist.sh ' + outdir + ' ' + prefix + ' ' + str(numcell) +
        ' ' + str(pv) + ' ' + str(FakeProb) + ' ' + gtffilepath + ' ' + codedir + '/../data/GenePos.txt ' + codedir)
    print('--------------------------------------------')
    if Rename:
        print('Please run the commands below to rename the files to their original names\n')
        print('python ' + codedir + 'RenameFastqFiles.py ' + filedir + ' ' + filedir + '/RenameList.txt')
except getopt.GetoptError:
    help()
