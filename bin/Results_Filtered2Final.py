from __future__ import print_function
from __future__ import division
import sys
import math
import os
import numpy

# ***** readme *****
# 这个代码从结果文件AllResults_filtered.txt中提取出概率低、p值低的出来
# 并且做出筛选

resultfile = open(sys.argv[1])
Pvaluethres = float(sys.argv[2])
probthres = float(sys.argv[3])
totalcellnum = min(int(sys.argv[4]), 500)


resultfilelines = resultfile.readlines()
goodpv = []
badpv = []
FDRDict = {}
setpcutoff = False
threspv = 20
setthrespv = False
for line in resultfilelines:
    if line.startswith('Fusion'):
        continue
    info = line.rstrip().split('\t')
    if (int(info[2]) / int(info[1])) < 1.25 or int(info[1]) < 5:
        badpv.append(float(info[7]))
    else:
        goodpv.append(float(info[7]))

for l in range(1000):
    pcutoff = pow(10, -l - 0.5)
    smallgoodpv = 0
    smallbadpv = 0
    for p in goodpv:
        if p < pcutoff:
            smallgoodpv += 1
    for p in badpv:
        if p < pcutoff:
            smallbadpv += 1
    aaa = smallbadpv / max(1, len(badpv)) * (len(goodpv) + len(badpv)) / max(1, smallgoodpv + smallbadpv)
    #sys.stderr.write(str(aaa) + '\n')
    if aaa < Pvaluethres and not setpcutoff:
        threspv = l + 0.5
        FDRDict[l] = aaa
        FDRDict[-1] = aaa
        setpcutoff = True
        continue
    if setpcutoff:
        FDRDict[l] = min(aaa, FDRDict[-1])
        FDRDict[-1] = FDRDict[l]
    


if len(sys.argv) > 5:
    threspv = float(sys.argv[5])
    for i in range(1000):
        FDRDict[i] = pow(10, -i)


for line in resultfilelines:
    if line.startswith('Fusion'):
        continue
    info = line.rstrip().split('\t')
    if float(info[6]) <= probthres and float(info[7]) <= pow(10, -threspv):
        if float(info[7]) > 0:
            fdrlevel = FDRDict[math.floor(-math.log10(float(info[7])) - 0.5)]
        else:
            fdrlevel = FDRDict[-1]
        info.insert(8, str(fdrlevel))
        sepe = '\t'
        print(sepe.join(info))
