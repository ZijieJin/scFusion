from __future__ import print_function
from __future__ import division
import sys
import math
import os
import numpy


# ***** readme *****
# This code finds the distribution of the chimeric read in each cell for each gene.


def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


def chr2num(str):
    if str.startswith('chr'):
        str = str[3:]
    if str == 'X':
        return 23
    if str == 'Y':
        return 24
    return int(str)


def SolveClip(str, readlen):
    num = []
    alp = []
    res = []
    start = 0
    scount = 0
    Ndel = 0
    length = len(str)
    lastalp = -1
    for i in range(length):
        if str[i].isalpha():
            alp.append(str[i])
            num.append(int(str[lastalp + 1:i]))
            lastalp = i
    i = 0
    while i < len(alp):
        if alp[i] == 'N':
            for j in range(i):
                if alp[j] == 'S':
                    scount = 1
                    break
            del alp[i]
            Ndel += num[i]
            del num[i]
            continue
        i += 1
    lastm = -1
    i = 0
    while i < len(alp):
        if alp[i] == 'M':
            if lastm == -1:
                lastm = i
                continue
            if lastm == i - 1:
                num[i - 1] += num[i]
                del alp[i]
                del num[i]
                continue
            lastm = i
        i += 1
    numsum = [num[0]]
    for i in range(1, len(num)):
        numsum.append(numsum[-1] + num[i])
    for i in range(len(alp)):
        if alp[i] == 'M':
            if i == 0:
                res = [numsum[0]]
            elif i == len(alp) - 1:
                res.append(numsum[i - 1])
            else:
                res.append(numsum[i - 1])
                res.append(numsum[i])
    for i in range(len(alp)):
        if alp[i] == 'M':
            if i > 0:
                start = numsum[i - 1]
            break
    return [res, alp, start, 0, Ndel, scount]


def AddCount(genename, q=1.0):
    genename = genename.split(';')
    for gene in genename:
        gene = gene.split('.')[0]
        if gene in ExprDic:
            if ExprDic[gene][1] == i:
                ExprDic[gene][0][-1] += q
            else:
                for k1 in range(ExprDic[gene][1]+1, i):
                    ExprDic[gene][0].append(0)
                ExprDic[gene][0].append(q)
                ExprDic[gene][1] = i
        else:
            ExprDic[gene] = [[], i]
            for k1 in range(start, i):
                ExprDic[gene][0].append(0)
            ExprDic[gene][0].append(q)


def AddCount2(genename, q=1.0):
    genename = genename.split('.')[0]
    if genename in ChimericExprDic:
        ChimericExprDic[genename] += q
    else:
        ChimericExprDic[genename] = q


argnum = len(sys.argv)
filedir = sys.argv[1]
start = int(sys.argv[2])
last = int(sys.argv[3])
ExprDir = sys.argv[4]
HomologyResultFile = open(sys.argv[5])
FusionMatrix = {}
EncompassingMatrix = {}
FusionPos = {}
Clusters = {}
ClusterSize = {}
ClusterCount = {}
ExprDic = {}
ChimericExprDic = {}
GeneDictionary = {}
HomologyDic = {}
good_dist = 1
goodcount = 0
prefix = ''
if argnum == 7:
    prefix = sys.argv[6]
if prefix == '.':
    prefix = ''
readfile = open(filedir + '../ChiDist/' + prefix + 'FusionRead.txt', 'w')
# for finding reads near the brkpnts
CandidateList = []
CandidateVector = {}


for line in HomologyResultFile.readlines():
    if line[0] == '#':
        continue
    info = line.split('\t')
    try:
        fusionscore = float(info[2])
        chr1 = info[3]
        pos1 = info[5]
        chr2 = info[4].replace(' ', '')
        pos2 = info[6]
        homoscore = info[11]
        gccontent = [info[12], info[13]]
        cellcount = info[8]
        read1 = info[16]
        read2 = info[17].replace('\n', '').replace('\r', '')
        HomologyDic[chr1 + '\t' + chr2 + '\t' + pos1 + '\t' + pos2] = [homoscore, gccontent, cellcount, read1, read2]
    except:
        pass
HomologyResultFile.close()
for i in range(int(start), int(last) + 1):
    try:
        infile = open(filedir + '/' + str(i) + '_FusionSupport.txt')
        exprfile = open(ExprDir + '/' + str(i) + '.rpkm.txt')
        samfile = open(filedir + '/' + str(i) + '_geneanno.sam')
    except:
        continue
    foundcount = 0
    for key in CandidateVector:
        if CandidateVector[key]:
            foundcount += 1
    sys.stderr.write('Starting: ' + str(i) + '\tCandidate Size: ' + str(len(CandidateVector)) + '\tFound Size: ' + str(foundcount) + '\n')
    lastgene = ''
    totalrc = 0
    for line in exprfile.readlines():
        info = line.split('\t')
        gene = info[0]
        rc = int(info[1])
        if gene == lastgene:
            totalrc += rc
            continue
        if lastgene != '':
            AddCount(lastgene, totalrc)
        lastgene = gene
        totalrc = rc
    AddCount(lastgene, totalrc)
    lastname = ''
    lines = []
    badreadname = ''
    allsamfilelines = samfile.readlines()
    readlength = 1000
    for line in allsamfilelines:
        if line[0] == '#':
            continue
        info = line.split('\t')
        thisname = info[1]
        if thisname == badreadname:
            continue
        if thisname == lastname:
            lines.append(line)
        else:
            if len(lines) == 2:
                info1 = lines[0].split('\t')
                info2 = lines[1].split('\t')
                if info1[0] != info2[0]:
                    gene1 = info1[0]
                    gene2 = info2[0]
                    AddCount(gene1)
                    AddCount(gene2)
                    AddCount2(gene1)
                    AddCount2(gene2)
            elif len(lines) == 3:
                info1 = lines[0].split('\t')
                info2 = lines[1].split('\t')
                info3 = lines[2].split('\t')
                gene1 = info1[0].split('.')[0]
                gene2 = info2[0].split('.')[0]
                gene3 = info3[0].split('.')[0]
                if gene1 != gene2 or gene2 != gene3:
                    readlength = len(info1[10])
                    if info1[10] == info2[10] or info1[10] == ReverseComplement(info2[10]):
                        AddCount(gene1, 0.5)
                        AddCount(gene2, 0.5)
                        AddCount(gene3)
                        AddCount2(gene1, 0.5)
                        AddCount2(gene2, 0.5)
                        AddCount2(gene3)
                        readinfo1 = info1
                        readinfo2 = info2
                    elif info1[10] == info3[10] or info1[10] == ReverseComplement(info3[10]):
                        AddCount(gene1, 0.5)
                        AddCount(gene2)
                        AddCount(gene3, 0.5)
                        AddCount2(gene1, 0.5)
                        AddCount2(gene2)
                        AddCount2(gene3, 0.5)
                        readinfo1 = info1
                        readinfo2 = info3
                    elif info2[10] == info3[10] or info2[10] == ReverseComplement(info3[10]):
                        AddCount(gene1)
                        AddCount(gene2, 0.5)
                        AddCount(gene3, 0.5)
                        AddCount2(gene1)
                        AddCount2(gene2, 0.5)
                        AddCount2(gene3, 0.5)
                        readinfo1 = info3
                        readinfo2 = info2
                    else:
                        sys.stderr.write('!!!!!' + info1[1])
                        lines = [line]
                        lastname = thisname
                        readinfo2 = [1]
                        readinfo1 = [1]  # don't run codes below, just want to read lines
            lines = [line]
            lastname = thisname

    CandidateList = []
    for line in infile.readlines():
        info = line.split('\t')
        gene1 = info[0]
        gene2 = info[1]
        chromo1 = info[4]
        chromo2 = info[5]
        encompass = int(info[2])
        if gene1 + '\t' + gene2 in EncompassingMatrix:
            if encompass > 0:
                EncompassingMatrix[gene1 + '\t' + gene2].append(encompass)
        elif gene2 + '\t' + gene1 in EncompassingMatrix:
            if encompass > 0:
                EncompassingMatrix[gene2 + '\t' + gene1].append(encompass)
        else:
            if encompass > 0:
                EncompassingMatrix[gene1 + '\t' + gene2] = [encompass]
            else:
                EncompassingMatrix[gene1 + '\t' + gene2] = []
        splitcount = int(info[3])
        if splitcount == 0:
            continue
        splitreadinfo = info[6].split(';')
        splitreadinfo = list(set(splitreadinfo))        # remove duplication
        subscore = []
        Pos = {}
        subgrp = []
        for item in splitreadinfo:
            if len(item) > 3:
                iteminfo = item.split('+')
                if not 6 <= int(iteminfo[1]) <= readlength - 6:
                    continue
                pos = iteminfo[0].split(',')
                if pos[0] + ',' + pos[1] in Pos:
                    Pos[pos[0] + ',' + pos[1]] += 1
                else:
                    Pos[pos[0] + ',' + pos[1]] = 1
                badmap = iteminfo[2].split(',')
                # editscore = max(50 - abs(50 - int(iteminfo[1])) - int(badmap[0]) - int(badmap[1]), 0)
                ssscore = 1
                if len(subgrp) == 0:
                    subgrp = [{pos[0] + ',' + pos[1]: 1}]
                    subscore.append(ssscore)
                else:
                    fin = False
                    for j in range(len(subgrp)):
                        for key in subgrp[j]:
                            grpos0 = key.split(',')[0]
                            grpos1 = key.split(',')[1]
                            if abs(int(pos[0]) - int(grpos0)) <= 3 * good_dist and abs(int(pos[1]) - int(grpos1)) <= 3 * good_dist:
                                fin = True
                                if pos[0] + ',' + pos[1] in subgrp[j]:
                                    subgrp[j][pos[0] + ',' + pos[1]] += 1
                                else:
                                    subgrp[j][pos[0] + ',' + pos[1]] = 1
                                subscore[j] += ssscore
                                break
                        if fin:
                            break
                    if not fin:
                        subgrp.append({pos[0] + ',' + pos[1]: 1})
                        subscore.append(ssscore)
        if len(subgrp) == 0:
            continue
        rawscore = []
        for j in subscore:
            if j > 0.5:
                rawscore.append(j)
            else:
                rawscore.append(0)
        rev = False
        AvePos = []
        for mdict in subgrp:  # calculate average pos for each cluster
            left = 0
            right = 0
            cc = 0
            for key in mdict:
                pos = key.split(',')
                left += int(pos[0]) * mdict[key]
                right += int(pos[1]) * mdict[key]
                cc += mdict[key]
            left = int(left / cc)
            right = int(right / cc)
            AvePos.append([left, right, cc])
        if gene1 + '\t' + gene2 in FusionMatrix:
            added = list(numpy.zeros(len(FusionMatrix[gene1 + '\t' + gene2])))
            for j in range(len(subgrp)):
                dist = []
                for k in FusionPos[gene1 + '\t' + gene2][1]:
                    d1 = abs(k[1][0] - AvePos[j][0])
                    d2 = abs(k[1][1] - AvePos[j][1])
                    if d1 <= good_dist and d2 <= good_dist:
                        dist.append(d1 + d2)
                    else:
                        dist.append(10 * good_dist)
                if min(dist) >= 3 * good_dist:
                    FusionPos[gene1 + '\t' + gene2][1].append([subgrp[j], AvePos[j], [i]])
                    FusionMatrix[gene1 + '\t' + gene2].append([rawscore[j]])
                    added.append(0)
                else:
                    smallindex = dist.index(min(dist))
                    if added[smallindex] == 0:
                        FusionMatrix[gene1 + '\t' + gene2][smallindex].append(rawscore[j])
                        added[smallindex] = 1
                    else:
                        FusionMatrix[gene1 + '\t' + gene2][smallindex][-1] += (rawscore[j])
                    old = FusionPos[gene1 + '\t' + gene2][1][smallindex][1]
                    newcc = AvePos[j][2] + old[2]
                    newleft = int((old[2] * old[0] + AvePos[j][2] * AvePos[j][0]) / newcc)
                    newright = int((old[2] * old[1] + AvePos[j][2] * AvePos[j][1]) / newcc)
                    FusionPos[gene1 + '\t' + gene2][1][smallindex][1] = [newleft, newright, newcc]
                    FusionPos[gene1 + '\t' + gene2][1][smallindex][2].append(i)
                if not [gene1, gene2, chromo1, chromo2, AvePos[j][0], AvePos[j][1]] in CandidateList and not [gene2, gene1, chromo2, chromo1, AvePos[j][1], AvePos[j][0]] in CandidateList:
                    CandidateList.append([gene1, gene2, chromo1, chromo2, AvePos[j][0], AvePos[j][1]])
            gene = gene1 + '\t' + gene2
        elif gene2 + '\t' + gene1 in FusionMatrix:
            added = list(numpy.zeros(len(FusionMatrix[gene2 + '\t' + gene1])))
            for j in range(len(subgrp)):
                dist = []
                for k in FusionPos[gene2 + '\t' + gene1][1]:
                    d1 = abs(k[1][0] - AvePos[j][1])
                    d2 = abs(k[1][1] - AvePos[j][0])
                    if d1 <= good_dist and d2 <= good_dist:
                        dist.append(d1 + d2)
                    else:
                        dist.append(10 * good_dist)
                if min(dist) >= 3 * good_dist:
                    FusionPos[gene2 + '\t' + gene1][1].append(
                        [subgrp[j], [AvePos[j][1], AvePos[j][0], AvePos[j][2]], [i]])
                    FusionMatrix[gene2 + '\t' + gene1].append([rawscore[j]])
                    added.append(0)
                else:
                    smallindex = dist.index(min(dist))
                    if added[smallindex] == 0:
                        FusionMatrix[gene2 + '\t' + gene1][smallindex].append(rawscore[j])
                        added[smallindex] = 1
                    else:
                        FusionMatrix[gene2 + '\t' + gene1][smallindex][-1] += (rawscore[j])
                    old = FusionPos[gene2 + '\t' + gene1][1][smallindex][1]
                    newcc = AvePos[j][2] + old[2]
                    newleft = int((old[2] * old[0] + AvePos[j][2] * AvePos[j][1]) / newcc)
                    newright = int((old[2] * old[1] + AvePos[j][2] * AvePos[j][0]) / newcc)
                    FusionPos[gene2 + '\t' + gene1][1][smallindex][1] = [newleft, newright, newcc]
                    FusionPos[gene2 + '\t' + gene1][1][smallindex][2].append(i)
                if not [gene1, gene2, chromo1, chromo2, AvePos[j][0], AvePos[j][1]] in CandidateList and not [gene2, gene1, chromo2, chromo1, AvePos[j][1], AvePos[j][0]] in CandidateList:
                    CandidateList.append([gene1, gene2, chromo1, chromo2, AvePos[j][0], AvePos[j][1]])
            gene = gene2 + '\t' + gene1
            rev = True
        else:
            FusionMatrix[gene1 + '\t' + gene2] = []
            for item in rawscore:
                FusionMatrix[gene1 + '\t' + gene2].append([item])
            FusionPos[gene1 + '\t' + gene2] = [[chromo1, chromo2], []]
            for k in range(len(subgrp)):
                FusionPos[gene1 + '\t' + gene2][1].append([subgrp[k], AvePos[k], [i]])
                if not [gene1, gene2, chromo1, chromo2, AvePos[k][0], AvePos[k][1]] in CandidateList and not [gene2, gene1, chromo2, chromo1, AvePos[k][1], AvePos[k][0]] in CandidateList:
                    CandidateList.append([gene1, gene2, chromo1, chromo2, AvePos[k][0], AvePos[k][1]])
            gene = gene1 + '\t' + gene2

    for l in range(len(CandidateList)):
        thischr1 = chr2num(CandidateList[l][2])
        thischr2 = chr2num(CandidateList[l][3])
        pos1 = int(CandidateList[l][4])
        pos2 = int(CandidateList[l][5])
        if thischr1 > thischr2 or (thischr1 == thischr2 and pos2 > pos1):
            thischr1, thischr2 = thischr2, thischr1
            pos1, pos2 = pos2, pos1
        if not (thischr1, thischr2, pos1, pos2) in CandidateVector:
            CandidateVector[(thischr1, thischr2, pos1, pos2)] = []

    for line in allsamfilelines:
        if line[0] == '#':
            continue
        info = line.split('\t')
        thisname = info[1]
        if thisname == badreadname:
            continue
        if thisname == lastname:
            lines.append(line)
        else:
            if len(lines) == 2:
                info1 = lines[0].split('\t')
                info2 = lines[1].split('\t')
            elif len(lines) == 3:
                info1 = lines[0].split('\t')
                info2 = lines[1].split('\t')
                info3 = lines[2].split('\t')
                gene1 = info1[0].split('.')[0]
                gene2 = info2[0].split('.')[0]
                gene3 = info3[0].split('.')[0]
                if gene1 != gene2 or gene2 != gene3:
                    if info1[10] == info2[10] or info1[10] == ReverseComplement(info2[10]):
                        readinfo1 = info1
                        readinfo2 = info2
                    elif info1[10] == info3[10] or info1[10] == ReverseComplement(info3[10]):
                        readinfo1 = info1
                        readinfo2 = info3
                    elif info2[10] == info3[10] or info2[10] == ReverseComplement(info3[10]):
                        readinfo1 = info3
                        readinfo2 = info2
                    else:
                        # sys.stderr.write('!!!!!' + info1[1])
                        lines = [line]
                        lastname = thisname
                        readinfo2 = [1]
                        readinfo1 = [1]  # don't run codes below, just want to read lines
                    if len(readinfo1) > 1 and len(gene1) > 1 and len(gene2) > 1 and len(gene3) > 1:  # normal conditions
                        gene1 = readinfo1[0]
                        gene2 = readinfo2[0]
                        chromo1 = readinfo1[3]
                        chromo2 = readinfo2[3]
                        if not chromo1.startswith('chr'):
                            chromo1 = 'chr' + chromo1
                        if not chromo2.startswith('chr'):
                            chromo2 = 'chr' + chromo2
                        clip1 = readinfo1[6]
                        clip2 = readinfo2[6]
                        readlen = len(readinfo1[10])
                        clipsplit1 = SolveClip(clip1, readlen)
                        clipsplit2 = SolveClip(clip2, readlen)
                        cc = False
                        splitpnt1 = -1
                        splitpnt2 = -1
                        for i in range(len(clipsplit1[0])):
                            for j in range(len(clipsplit2[0])):
                                if clipsplit2[0][j] == clipsplit1[0][i] or clipsplit2[0][j] + clipsplit1[0][i] == readlen:
                                    splitpnt1 = clipsplit1[0][i]
                                    splitpnt2 = clipsplit2[0][j]
                                    cc = True
                                    break
                            if cc:
                                break
                        if readlen - 0.5 > splitpnt1 > -0.5 and readlen - 0.5 > splitpnt2 > -0.5:
                            if clipsplit1[1][i] == 'S':
                                brkpnt1 = int(readinfo1[4])
                                direct1 = '+'
                            else:
                                brkpnt1 = splitpnt1 - clipsplit1[2] + int(readinfo1[4]) + clipsplit1[4] - 1
                                direct1 = '-'
                            if clipsplit2[1][j] == 'S':
                                brkpnt2 = int(readinfo2[4])
                                direct2 = '+'
                            else:
                                brkpnt2 = splitpnt2 - clipsplit2[2] + int(readinfo2[4]) + clipsplit2[4] - 1
                                direct2 = '-'
                            thisposlist = [chr2num(chromo1), chr2num(chromo2), brkpnt1, brkpnt2]
                            exchange = False
                            if thisposlist[0] > thisposlist[1] or (thisposlist[0] == thisposlist[1] and thisposlist[3] > thisposlist[2]):
                                thisposlist = [thisposlist[1], thisposlist[0], thisposlist[3], thisposlist[2]]
                                exchange = True
                            found = False
                            for k1 in range(-3, 4):
                                for k2 in range(-3, 4):
                                    if (thisposlist[0], thisposlist[1], thisposlist[2] + k1, thisposlist[3] + k2) in CandidateVector:
                                        found = True
                                        thisposindex = (thisposlist[0], thisposlist[1], thisposlist[2] + k1, thisposlist[3] + k2)
                            if found:
                                extractedRead = ''
                                if CandidateVector[thisposindex]:
                                    if CandidateVector[thisposindex][1] != 30:
                                        dist = min(CandidateVector[thisposindex][1], 60 - CandidateVector[thisposindex][1])
                                        if 30 > splitpnt1 > dist:
                                            extractedRead = readinfo1[10][:60]
                                            inbrkpnt = splitpnt1
                                        if 30 <= splitpnt1 <= readlen - 30:
                                            extractedRead = readinfo1[10][splitpnt1 - 30:splitpnt1 + 30]
                                            inbrkpnt = 30
                                        if readlen - dist > splitpnt1 > readlen - 30:
                                            extractedRead = readinfo1[10][readlen - 60:]
                                            inbrkpnt = splitpnt1 - readlen + 60
                                        if extractedRead != '':
                                            if direct1 == '+':
                                                extractedRead = ReverseComplement(extractedRead)
                                                inbrkpnt = 60 - inbrkpnt
                                            if exchange:
                                                extractedRead = ReverseComplement(extractedRead)
                                                inbrkpnt = 60 - inbrkpnt
                                            CandidateVector[thisposindex][0] = extractedRead
                                            CandidateVector[thisposindex][1] = inbrkpnt
                                else:
                                    if splitpnt1 < 30:
                                        extractedRead = readinfo1[10][:60]
                                        inbrkpnt = splitpnt1
                                    if 30 <= splitpnt1 <= readlen - 30:
                                        extractedRead = readinfo1[10][splitpnt1 - 30:splitpnt1 + 30]
                                        inbrkpnt = 30
                                    if splitpnt1 > readlen - 30:
                                        extractedRead = readinfo1[10][readlen - 60:]
                                        inbrkpnt = splitpnt1 - readlen + 60
                                    if direct1 == '+':
                                        extractedRead = ReverseComplement(extractedRead)
                                        inbrkpnt = 60 - inbrkpnt
                                    if exchange:
                                        extractedRead = ReverseComplement(extractedRead)
                                        inbrkpnt = 60 - inbrkpnt
                                        direct1, direct2 = direct2, direct1
                                    CandidateVector[thisposindex] = [extractedRead, inbrkpnt, direct1, direct2]
            lines = [line]
            lastname = thisname
    samfile.close()
    infile.close()
    exprfile.close()
for genename in ExprDic:
    if ExprDic[genename][1] < last:
        for k in range(ExprDic[genename][1]+1, last+1):
            ExprDic[genename][0].append(0)
for gene in FusionMatrix:
    temp = FusionPos[gene][1]
    for i in range(len(temp)):
        keylist = sorted(temp[i][0], key=temp[i][0].__getitem__, reverse=True)
        possum = 0
        countlist = []
        for key in keylist:
            possum += temp[i][0][key]
            countlist.append(int(temp[i][0][key]))
        if possum > 0 and max(countlist) / possum >= 0.4:
            for key in temp[i][0]:
                if key == keylist[0]:
                    pos0 = int(keylist[0].split(',')[0])
                    pos1 = int(keylist[0].split(',')[1])
                    if abs(pos0 - temp[i][1][0]) + abs(pos1 - temp[i][1][1]) < abs(pos1 - temp[i][1][0]) + abs(
                            pos0 - temp[i][1][1]):
                        FusionPos[gene][1][i].append([pos0, pos1])
                    else:
                        FusionPos[gene][1][i].append([pos1, pos0])
                    break
        else:
            FusionPos[gene][1][i].append([-1, -1])
templines = []
for gene in FusionMatrix:
    for i in range(len(FusionPos[gene][1])):
        if FusionPos[gene][1][i][3][0] > -0.5:
            genepart = gene.split('\t')
            usecell = FusionPos[gene][1][i][2]
            usecell = list(set(usecell))
            try:
                aaaaa = []
                thisgenename = genepart[0].split('.')[0]
                for cell in usecell:
                    aaaaa.append(ExprDic[thisgenename][0][cell-start])
            except:
                aaaaa = list(numpy.zeros(last-start+1, int))
                continue
            try:
                bbbbb = []
                thisgenename = genepart[1].split('.')[0]
                for cell in usecell:
                    bbbbb.append(ExprDic[thisgenename][0][cell-start])
            except:
                bbbbb = list(numpy.zeros(last-start+1, int))
                continue
            if FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1] + '\t' + str(FusionPos[gene][1][i][3][0]) + '\t' \
                    + str(FusionPos[gene][1][i][3][1]) in HomologyDic:
                [homoscore, gccontent, cellcount, read1, read2] = HomologyDic[FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1]
                                                                + '\t' + str(FusionPos[gene][1][i][3][0]) + '\t' + str(
                    FusionPos[gene][1][i][3][1])]
            elif FusionPos[gene][0][1] + '\t' + FusionPos[gene][0][0] + '\t' + str(FusionPos[gene][1][i][3][1]) + '\t' \
                    + str(FusionPos[gene][1][i][3][0]) in HomologyDic:
                [homoscore, gccontent, cellcount, read1, read2] = HomologyDic[FusionPos[gene][0][1] + '\t' + FusionPos[gene][0][0]
                                                                + '\t' + str(FusionPos[gene][1][i][3][1]) + '\t' + str(
                    FusionPos[gene][1][i][3][0])]
            else:
                continue
            engene = gene
            if gene not in EncompassingMatrix:
                engene = genepart[1] + '\t' + genepart[0]
            templine = genepart[0] + '\t' + genepart[1] + '\t' + \
                       str(len(FusionMatrix[gene][i])) + '\t' + str(FusionMatrix[gene][i]) + '\t' + \
                       str(len(EncompassingMatrix[engene])) + '\t' + str(EncompassingMatrix[engene]) + '\t' + \
                       FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1] + '\t' + str(FusionPos[gene][1][i][3][0]) + \
                       '\t' + str(FusionPos[gene][1][i][3][1]) + '\t' + str(i) + '\t' + homoscore + '\t' + \
                       gccontent[0] + '\t' + gccontent[1] + '\t' + str(aaaaa) + '\t' + str(bbbbb) + '\t' + read1 + '\t' + read2
            templines.append(templine)
uselines = []
for i in range(len(templines)):
    info = templines[i].split('\t')
    thischr1 = chr2num(info[6])
    thischr2 = chr2num(info[7])
    thispos1 = int(info[8])
    thispos2 = int(info[9])
    exchange = False
    if thischr1 > thischr2 or (thischr1 == thischr2 and thispos2 > thispos1):
        thischr1, thischr2 = thischr2, thischr1
        thispos1, thispos2 = thispos2, thispos1
        exchange = True
    found = False
    for k1 in range(-3, 4):
        for k2 in range(-3, 4):
            if (thischr1, thischr2, thispos1 + k1, thispos2 + k2) in CandidateVector:
                found = True
                thisposindex = (thischr1, thischr2, thispos1 + k1, thispos2 + k2)
    if found:
        if CandidateVector[thisposindex]:
            if exchange:
                uselines.append(templines[i].rstrip() + '\t' + ReverseComplement(CandidateVector[thisposindex][
                    0]) + '\t' + str(60-CandidateVector[thisposindex][1]) + '\t' + CandidateVector[thisposindex][3] + '\t' +
                                CandidateVector[thisposindex][2])
            else:
                uselines.append(templines[i].rstrip() + '\t' + CandidateVector[thisposindex][
                    0] + '\t' + str(CandidateVector[thisposindex][1]) + '\t' + CandidateVector[thisposindex][2] + '\t' +
                                CandidateVector[thisposindex][3])


linedic = {}
for line in uselines:
    info = line.split('\t')
    score = len(info[3])
    linedic[line] = score
for key in sorted(linedic, key=linedic.__getitem__, reverse=True):
    info = key.split('\t')
    read = info[-4]
    splitpos = info[-3]
    chromo1 = info[6]
    chromo2 = info[7]
    pos1 = info[8]
    pos2 = info[9]
    if read.find('N') > -1:
        continue
    if len(read) == 60 or True:
        readfile.write(read + '\t' + splitpos + '\t' + chromo1 + ':' + pos1 + ':' + info[-2] + '\t' + chromo2 + ':' +
                       pos2 + ':' + info[-1] + '\n')
        print(key)
readfile.close()

