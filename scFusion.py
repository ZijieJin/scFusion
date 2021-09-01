from __future__ import print_function
from __future__ import division
import sys
import getopt
import os
import threading
import numpy
import subprocess


def help():
    print('----------Command Line Options: ----------')
    print('Parameters must be specified: ')
    print('-f, --FileDir: The folder of single cell sequencing files')
    print('-b, --Begin: The first index of single cell sequencing file')
    print('-e, --End: The last index of single cell sequencing file')
    print('-s, --STARReference: The reference folder of STAR')
    print('Parameters with default value, but should be changed for your setting: ')
    print('-g, --Genome: The genome reference file (*.fasta or *.fa), default is hg19.fa in the data folder')
    print('-a, --Annotation: The gtf annotation file (*.gtf), default is ref_annot.gtf in the data folder')
    print('-m, --Mappability: The mappability file, default is hg19mappability75.txt in the data folder, left blank for keeping all reads')
    print('-o, --OutDir: The output folder of the results, default is the same as FileDir')
    print('-t, --Thread: Number of threads, at least 4, default is 8')
    print('-l, --LimitThread: Number of maximum threads allowed for each STAR mapping task, default is 20')
    print('-w, --Weight: The weight file of the deep-learning network, default is weight-V9-2.hdf5 in the data folder')
    print('-E, --Epoch: The number of epochs in the retraining step between 3 and 999, default is 10')
    print('-p, --Prefix: The prefix of result file, default is blank')
    print('-v, --PvalueCutoff: Pvalue cutoff, default is 0.05')
    print('-n, --NetworkCutoff: Network score cutoff, default is 0.75')
    print('Step Controls:')
    print('--Rename: Rename the files with the consecutive indexes')
    print('--SkipMapping: Skip STAR Mapping, if you already have the mapping result at OutDir/StarMapping/')
    print('--SkipBS: Skip the basic processing step')
    print('--SkipCombining: Skip the combining step')
    print('--SkipRetrain: Skip the retraining step, and apply the weights specified by -w to the network')
    print('--SkipPredict: Skip the predicting step using network')


def printlog(str):
    print(str)
    logfile.write(str)



def RunSTARMapping(s, e, core):
    os.system('sh ' + codedir + 'StarMapping_Chimeric.sh ' + filedir + ' ' + str(s) + ' ' + str(e) + ' ' + outdir +
              'StarMapping/ ' + starrefdir + ' ' + str(core))
    printlog('Finish mapping! Index: ' + str(s) + ' ~ ' + str(e) + '\n')


def RunBS(s, e):
    os.system('sh ' + codedir + 'CombinePipeline_before_FS.sh ' + outdir + ' ' + str(s) + ' ' + str(e) + ' ' +
              gtffilepath + ' ' + mappabilityfilepath + ' ' + exonposfilepath + ' ' + codedir)
    printlog('Finish Basic Processing! Index: ' + str(s) + ' ~ ' + str(e) + '\n')


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
    epoch = 10
    prefix = '.'
    start = -10000000
    end = -10000000
    cellindex = []
    opts, args = getopt.getopt(sys.argv[1:], 'ht:s:g:a:b:e:f:m:p:w:v:n:o:l:E:',
                               ['help', 'Thread=', 'STARReference=', 'Genome=', 'Annotation=', 'Begin=', 'End=',
                                'FileDir=', 'Mappability=', 'LimitThread=',
                                'Prefix=', 'Weight=', 'Epoch=', 'PvalueCutoff=', 'NetworkCutoff=', 'OutDir=', 'SkipMapping',
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
        elif opt in ('-E', '--Epoch'):
            epoch = int(arg)
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
    logfile = open(outdir + 'log.txt', 'w')
    printlog('\nGenerating some necessary files!\n')
    if not os.path.exists(gtffilepath + '.added'):
        os.system('python ' + codedir + 'Addchr2gtf.py ' + gtffilepath + ' > ' + gtffilepath + '.added')
    gtffilepath = gtffilepath + '.added'
    if not os.path.exists(codedir + '/../data/GenePos.txt'):
        os.system('python ' + codedir + 'GetGenePos.py ' + gtffilepath + ' > ' + codedir + '/../data/GenePos.txt')
    if not os.path.exists(exonposfilepath):
        os.system('python ' + codedir + 'GetExonPos.py ' + gtffilepath + ' > ' + exonposfilepath)
    for i in range(start, end + 1):
        if os.path.exists(filedir + str(i) + '_2.fastq'):
            cellindex.append(i)
    if not SkipBS and not os.path.exists(gtffilepath[:-5] + 'db'):
        aaa = subprocess.check_output(
                'pyensembl install --reference-name GRCH37 --annotation-name my_genome_features --gtf ' + gtffilepath,
                shell=True, stderr=subprocess.STDOUT)
        logfile.write(str(aaa))
    numcell = len(cellindex)
    printlog('Parameter Check Complete!\n')
    os.system('mkdir -p ' + outdir + '/StarMapping/')
    os.system('mkdir -p ' + outdir + '/ChimericOut/')
    os.system('mkdir -p ' + outdir + '/Expr/')
    printlog(str(numcell) + ' cell sequencing files found!\n')
    if Rename:
        printlog('Start to rename files!')
        aaa = subprocess.check_output(
            'python ' + codedir + 'RenameFastqFiles.py ' + filedir, stderr=subprocess.STDOUT, shell=True)
        logfile.write(str(aaa))
        printlog('Successfully rename files to consecutive indexes!')
    if not SkipMapping:
        if numthread <= 24:
            printlog(
                'Start mapping! Index: ' + str(start) + ' ~ ' + str(end) + ', using core: ' + str(numthread - 1) + '\n')
            RunSTARMapping(start, end, numthread - 1)
        else:
            numtask = min(round(numthread / min(16, maxthread)), numcell)
            numcore = int(min(numpy.floor((numthread - 1) / numtask), maxthread))
            numcelleachtask = int(numpy.ceil(numcell / numtask))
            threads = []
            for j in range(numtask):
                printlog('Start mapping! Index: ' + str(cellindex[j * numcelleachtask]) + ' ~ ' + str(
                    min(cellindex[(j + 1) * numcelleachtask - 1], end)) + ', using core: ' + str(numcore) + '\n')
                threads.append(threading.Thread(target=RunSTARMapping, args=(
                cellindex[j * numcelleachtask], min(cellindex[min((j + 1) * numcelleachtask - 1, numcell - 1)], end), numcore)))
            for t in threads:
                t.start()
            while True:
                os.system('sleep 1')
                flag = True
                for t in threads:
                    if t.is_alive():
                        flag = False
                if flag:
                    break
    else:
        printlog('Mapping Skipped\n')
    if not SkipBS:
        numtask = min(numthread - 1, numcell)
        numcelleachtask = int(numpy.ceil(numcell / numtask))
        actualnumtask = int(numpy.ceil(numcell / numcelleachtask))
        threads = []
        for j in range(actualnumtask):
            printlog('Start Basic Processing! Index: ' + str(cellindex[j * numcelleachtask]) + ' ~ ' + str(
                min(cellindex[(j + 1) * numcelleachtask - 1], end)) + ', using core: 1\n')
            threads.append(threading.Thread(target=RunBS, args=(
            cellindex[j * numcelleachtask], min(cellindex[(j + 1) * numcelleachtask - 1], end))))
        for t in threads:
            t.start()
        while True:
            os.system('sleep 1')
            flag = True
            for t in threads:
                if t.is_alive():
                    flag = False
            if flag:
                break
    else:
        printlog('Basic Processing Skipped\n')
    if not SkipCombining:
        printlog('Start Combining!\n')
        aaa = subprocess.check_output(
            'sh ' + codedir + 'CombinePipeline_startwith_FS.sh ' + outdir + ' ' + str(start) + ' ' + str(
                end) + ' ' + prefix + ' ' +
            weightfilepath + ' ' + hg19filepath + ' ' + gtffilepath + ' ' + codedir, stderr=subprocess.STDOUT,
            shell=True)
        logfile.write(str(aaa))
        printlog('Finish Combining!\n')
    else:
        printlog('Combining Skipped\n')
    if not SkipRetrain:
        printlog('Start Retraining the Neural Network!\n')
        aaa = subprocess.check_output(
            'sh ' + codedir + 'CombinePipeline_Retrain.sh ' + outdir + ' ' + prefix + ' ' + weightfilepath + ' ' + str(epoch) + ' ' + codedir, stderr=subprocess.STDOUT,
            shell=True)
        logfile.write(str(aaa))
        maxep = -1
        for file in os.listdir(outdir + '/weights/'):
            epnumber = int(file[14:17])
            if epnumber > maxep:
                maxep = epnumber
        if maxep < 0:
            sys.stderr.write('Cannot create retraining weight files\n')
            sys.exit()
        if maxep < 10:
            epnumberstr = '00' + str(maxep)
        elif maxep < 100:
            epnumberstr = '0' + str(maxep)
        else:
            epnumberstr = str(maxep)
        weightfilepath = outdir + '/weights/RetrainWeight-' + epnumberstr + '.hdf5'
        printlog('Finish Retraining!\n')
        printlog('Using ' + weightfilepath + ' as weight file!\n')
    else:
        printlog('Retraining Step Skipped\n')
    if not SkipPredict:
        printlog('Start Predicting!\n')
        aaa = subprocess.check_output(
            'sh ' + codedir + 'CombinePipeline_Predict.sh ' + outdir + ' ' + str(start) + ' ' + str(
                end) + ' ' + prefix + ' ' +
            weightfilepath + ' ' + hg19filepath + ' ' + gtffilepath + ' ' + codedir, stderr=subprocess.STDOUT,
            shell=True)
        logfile.write(str(aaa))
        printlog('Finish Predicting Step using Neural Network!\n')
    else:
        printlog('Predict Step Skipped\n')
    printlog('Start Statistical Model!\n')
    aaa = subprocess.check_output(
        'sh ' + codedir + 'CombinePipeline_startwith_ChiDist.sh ' + outdir + ' ' + prefix + ' ' + str(numcell) +
        ' ' + str(pv) + ' ' + str(FakeProb) + ' ' + gtffilepath + ' ' + codedir + '/../data/GenePos.txt ' + codedir,
        stderr=subprocess.STDOUT, shell=True)
    logfile.write(str(aaa))
    if Rename:
        printlog('Start to rename files!')
        aaa = subprocess.check_output(
            'python ' + codedir + 'RenameFastqFiles.py ' + filedir + ' ' + filedir + '/RenameList.txt', stderr=subprocess.STDOUT, shell=True)
        logfile.write(str(aaa))
        printlog('Successfully rename files to original names!')
    printlog('All Finished!\n')
except getopt.GetoptError:
    help()
