from __future__ import print_function
from __future__ import division
import sys
import numpy

# ***** readme *****
# This code filters the fusion candidates than have lower than 3 supportors.


def MakeString(mylist, sep='\t'):
    string = ''
    for i in range(len(mylist) - 1):
        string += str(mylist[i]) + sep
    string += mylist[-1]
    return string


bigcount = 0
ChiDistFile = open(sys.argv[1])
templines = []
linedic = {}
for line in ChiDistFile.readlines():
    if line[0] == '#':
        continue
    try:
        info = line.split('\t')
        cc = info[3].replace('\n', '').replace('\r', '')
        count = cc[1:-1].split(', ')
        for i in range(len(count) - 1, -1, -1):
            if int(count[i].replace("'", '')) < 1:
                count.remove(count[i].replace("'", ''))
        mysum = 0
        for item in count:
            mysum += int(item.replace("'", ''))
        if len(count) == 1:
            continue
        if mysum / len(count) > 5:
            continue
        info[3] = '[' + MakeString(count, ', ') + ']'
        info[2] = len(count)
        templines.append(MakeString(info))
    except:
        pass
ChiDistFile.close()
for line in templines:
    info = line.split('\t')
    score = len(info[3])
    linedic[line] = score
numrecord = len(linedic)
for key in sorted(linedic, key=linedic.__getitem__, reverse=True):
    if bigcount < min(10.0, numpy.floor(numrecord / 50)):
        bigcount += 1
    else:
        print(key, end='')
