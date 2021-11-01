from __future__ import print_function
from __future__ import division
import sys
import os.path


# ***** readme *****
# 输入ref_annot文件，输出一个列表，是基因的位置

reffile = open(sys.argv[1])
genenameset = []
for line in reffile.readlines():
    if line.startswith('#'):
        continue
    info = line.split('\t')
    if info[2] == 'gene':
        chromo = info[0]
        start = info[3]
        end = info[4]
        geneinfo = info[8].split('; ')
        genename = ''
        for item in geneinfo:
            if item.find('gene_name') > -1:
                genename = item[11:-1]
        if genename == '':
            for item in geneinfo:
                if item.find('gene_id') > -1:
                    genename = item[9:-1]
        if genename == '':
            continue
        if genename not in genenameset:
            print(chromo + '\t' + start + '\t' + end + '\t' + genename)
