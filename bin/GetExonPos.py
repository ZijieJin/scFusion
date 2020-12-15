from __future__ import print_function
from __future__ import division
import sys
import os.path


# ***** readme *****
# 输入ref_annot文件，输出一个列表，是基因的位置

reffile = open(sys.argv[1])
for line in reffile.readlines():
    if line.startswith('#'):
        continue
    info = line.split('\t')
    if info[2] == 'exon':
        chromo = info[0]
        start = info[3]
        end = info[4]
        geneinfo = info[8].split('; ')
        genename = ''
        genetype = ''
        for item in geneinfo:
            if item.find('gene_name') > -1:
                genename = item[11:-1]
            if item.find('gene_type') > -1:
                genetype = item[11:-1]
        if genename != '' and genetype != '':
            print(chromo + '\t' + start + '\t' + end + '\t' + genename)