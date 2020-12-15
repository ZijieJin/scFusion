from __future__ import print_function
from __future__ import division
import sys
import math
import os


# ***** readme *****
# This code calculate fusion score for each cell and
# add rawscore for every cell


def sigmoid(x):
    return 1 / (1 + math.exp(-x))


def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


def AddCount(genenum, q=1.0):
    genenum = genenum.split('.')[0]
    if genenum in ExprDic:
        ExprDic[genenum] += q
    else:
        ExprDic[genenum] = q


def GetExpr(gene):
    genelist = gene.split(',')
    sum = 0
    count = 0
    for item in genelist:
        if item in ExprDic:
            sum += ExprDic[item]
            count += 1
    if count > 0:
        return sum / count
    return 0


argnum = len(sys.argv)
filedir = sys.argv[1]
start = sys.argv[2]
last = sys.argv[3]
ExprDir = sys.argv[4]
FusionMatrix = {}
FusionPos = {}
Clusters = {}
ClusterSize = {}
ClusterCount = {}
ExprDic = {}
GeneDictionary = {}
good_dist = 20
EnableClustering = False
EmcmpscoreDic = {}
if argnum == 6:
    EnableClustering = True
    ClusteringFile = open(sys.argv[5])
    for line in ClusteringFile.readlines():
        if line[0] == '#':
            continue
        info = line.split('\t')
        cellname = info[0]
        cluster = info[1]
        Clusters[cellname] = cluster
        if cluster in ClusterSize:
            ClusterSize[cluster] += 1
        else:
            ClusterSize[cluster] = 1
        ClusterCount[cluster] = 0  # this var is for the final output fusion categories.
    ClusteringFile.close()
for i in range(int(start), int(last) + 1):
    try:
        infile = open(filedir + '/' + str(i) + '_FusionSupport.txt')
        exprfile = open(ExprDir + '/' + str(i) + '.rpkm.txt')
        samfile = open(filedir + '/' + str(i) + '_geneanno.sam')
    except:
        continue
    lastgene = ''
    totalrc = 0
    cellsup = 0
    ExprDic = {}
    for line in exprfile.readlines():
        info = line.split('\t')
        gene = info[0]
        rc = int(info[1])
        if gene == lastgene:
            totalrc += rc
            continue
        if lastgene != '':
            ExprDic[lastgene] = totalrc
        lastgene = gene
        totalrc = rc
    exprfile.close()
    ExprDic[lastgene] = totalrc
    lastname = ''
    lines = []
    badreadname = ''
    for line in samfile.readlines():
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
            elif len(lines) == 3:
                info1 = lines[0].split('\t')
                info2 = lines[1].split('\t')
                info3 = lines[2].split('\t')
                gene1 = info1[0].split('.')[0]
                gene2 = info2[0].split('.')[0]
                gene3 = info3[0].split('.')[0]
                if gene1 != gene2 or gene2 != gene3:
                    if info1[10] == info2[10] or info1[10] == ReverseComplement(info2[10]):
                        AddCount(gene1, 0.5)
                        AddCount(gene2, 0.5)
                        AddCount(gene3)
                    elif info1[10] == info3[10] or info1[10] == ReverseComplement(info3[10]):
                        AddCount(gene1, 0.5)
                        AddCount(gene2)
                        AddCount(gene3, 0.5)
                    elif info2[10] == info3[10] or info2[10] == ReverseComplement(info3[10]):
                        AddCount(gene1)
                        AddCount(gene2, 0.5)
                        AddCount(gene3, 0.5)
                    else:
                        sys.stderr.write('!!!!!' + info1[1])
                        lines = [line]
                        lastname = thisname
                        readinfo2 = [1]
                        readinfo1 = [1]  # don't run codes below, just want to read lines
            lines = [line]
            lastname = thisname
    samfile.close()
    for line in infile.readlines():
        info = line.split('\t')
        gene1 = info[0]
        gene2 = info[1]
        chromo1 = info[4]
        chromo2 = info[5]
        encompass = int(info[2])
        splitcount = int(info[3])
        if splitcount + encompass == 0:
            continue
        splitreadinfo = info[6].split(';')
        splitreadinfo = list(set(splitreadinfo))        # remove duplication
        encmpscore = encompass
        subscore = []
        Pos = {}
        try:
            scale1 = ExprDic[gene1]
            scale2 = ExprDic[gene2]
        except:
            sys.stderr.write('!!!' + gene1 + '\t' + gene2)
            continue
        scale = math.sqrt(math.log(1+scale1, 2) * math.log(1+scale2, 2))
        subgrp = []
        for item in splitreadinfo:
            if len(item) > 3:
                iteminfo = item.split('+')
                pos = iteminfo[0].split(',')
                if pos[0] + ',' + pos[1] in Pos:
                    Pos[pos[0] + ',' + pos[1]] += 1
                else:
                    Pos[pos[0] + ',' + pos[1]] = 1
                badmap = iteminfo[2].split(',')
                editscore = max(50 - abs(50 - int(iteminfo[1])) - int(badmap[0]) - int(badmap[1]), 0)
                if editscore > 12:
                    ssscore = 12 + (editscore - 12) / 10
                else:
                    ssscore = editscore
                if len(subgrp) == 0:
                    subgrp = [{pos[0] + ',' + pos[1]: 1}]
                    subscore.append(ssscore)
                else:
                    fin = False
                    for j in range(len(subgrp)):
                        for key in subgrp[j]:
                            grpos0 = key.split(',')[0]
                            grpos1 = key.split(',')[1]
                            if abs(int(pos[0]) - int(grpos0)) <= 20 and abs(int(pos[1]) - int(grpos1)) <= 20:
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
        rawscore = []
        for j in subscore:
            if j > 0.5 or encmpscore > 0:
                # rawscore.append(sigmoid((j + encmpscore) / 20 - 1.7))
                rawscore.append(pow(1 + (j + encmpscore) / (7 + scale), 1.3))
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
            left = left / cc
            right = right / cc
            AvePos.append([left, right, cc])
        if gene1 + '\t' + gene2 in FusionMatrix:
            if len(subgrp) == 0:
                if gene1 + '\t' + gene2 in EmcmpscoreDic:
                    EmcmpscoreDic[gene1 + '\t' + gene2] += rawscore
                elif gene2 + '\t' + gene1 in EmcmpscoreDic:
                    EmcmpscoreDic[gene2 + '\t' + gene1] += rawscore
                else:
                    EmcmpscoreDic[gene1 + '\t' + gene2] = rawscore
            for j in range(len(subgrp)):
                if len(FusionPos[gene1 + '\t' + gene2][1]) == 0 and len(subgrp) > 0:
                    FusionMatrix[gene1 + '\t' + gene2] = rawscore
                    for k in range(len(subgrp)):
                        FusionPos[gene1 + '\t' + gene2][1].append([subgrp[k], AvePos[k], [i], 1])
                else:
                    for subpos in AvePos:
                        dist = []
                        for k in FusionPos[gene1 + '\t' + gene2][1]:
                            d1 = abs(k[1][0] - subpos[0])
                            d2 = abs(k[1][1] - subpos[1])
                            if d1 <= good_dist and d2 <= good_dist:
                                dist.append(d1 + d2)
                            else:
                                dist.append(10 * good_dist)
                        if len(dist) > 0:
                            if min(dist) >= 3 * good_dist:
                                FusionPos[gene1 + '\t' + gene2][1].append([subgrp[j], AvePos[j], [i], 1])
                                FusionMatrix[gene1 + '\t' + gene2].append(rawscore[j])
                            else:
                                smallindex = dist.index(min(dist))
                                FusionMatrix[gene1 + '\t' + gene2][smallindex] += rawscore[j]
                                old = FusionPos[gene1 + '\t' + gene2][1][smallindex][1]
                                newcc = subpos[2] + old[2]
                                newleft = (old[2] * old[0] + subpos[2] * subpos[0]) / newcc
                                newright = (old[2] * old[1] + subpos[2] * subpos[1]) / newcc
                                FusionPos[gene1 + '\t' + gene2][1][smallindex][1] = [newleft, newright, newcc]
                                FusionPos[gene1 + '\t' + gene2][1][smallindex][2].append(i)
                                FusionPos[gene1 + '\t' + gene2][1][smallindex][3] += 1
            gene = gene1 + '\t' + gene2
        elif gene2 + '\t' + gene1 in FusionMatrix:
            if len(subgrp) == 0:
                if gene1 + '\t' + gene2 in EmcmpscoreDic:
                    EmcmpscoreDic[gene1 + '\t' + gene2] += rawscore
                elif gene2 + '\t' + gene1 in EmcmpscoreDic:
                    EmcmpscoreDic[gene2 + '\t' + gene1] += rawscore
                else:
                    EmcmpscoreDic[gene1 + '\t' + gene2] = rawscore
            for j in range(len(subgrp)):
                if len(FusionPos[gene2 + '\t' + gene1][1]) == 0 and len(subgrp) > 0:
                    FusionMatrix[gene2 + '\t' + gene1] = rawscore
                    for k in range(len(subgrp)):
                        FusionPos[gene2 + '\t' + gene1][1].append(
                            [subgrp[j], [AvePos[j][1], AvePos[j][0], AvePos[j][2]], [i], 1])
                else:
                    for subpos in AvePos:
                        dist = []
                        for k in FusionPos[gene2 + '\t' + gene1][1]:
                            d1 = abs(k[1][0] - subpos[1])
                            d2 = abs(k[1][1] - subpos[0])
                            if d1 <= good_dist and d2 <= good_dist:
                                dist.append(d1 + d2)
                            else:
                                dist.append(10 * good_dist)
                        if min(dist) >= 3 * good_dist:
                            FusionPos[gene2 + '\t' + gene1][1].append(
                                [subgrp[j], [AvePos[j][1], AvePos[j][0], AvePos[j][2]], [i], 1])
                            FusionMatrix[gene2 + '\t' + gene1].append(rawscore[j])
                        else:
                            smallindex = dist.index(min(dist))
                            FusionMatrix[gene2 + '\t' + gene1][smallindex] += rawscore[j]
                            old = FusionPos[gene2 + '\t' + gene1][1][smallindex][1]
                            newcc = subpos[2] + old[2]
                            newleft = (old[2] * old[0] + subpos[2] * subpos[1]) / newcc
                            newright = (old[2] * old[1] + subpos[2] * subpos[0]) / newcc
                            FusionPos[gene2 + '\t' + gene1][1][smallindex][1] = [newleft, newright, newcc]
                            FusionPos[gene2 + '\t' + gene1][1][smallindex][2].append(i)
                            FusionPos[gene2 + '\t' + gene1][1][smallindex][3] += 1
            gene = gene2 + '\t' + gene1
            rev = True
        else:
            if len(subgrp) == 0:
                if gene1 + '\t' + gene2 in EmcmpscoreDic:
                    EmcmpscoreDic[gene1 + '\t' + gene2] += rawscore
                elif gene2 + '\t' + gene1 in EmcmpscoreDic:
                    EmcmpscoreDic[gene2 + '\t' + gene1] += rawscore
                else:
                    EmcmpscoreDic[gene1 + '\t' + gene2] = rawscore
            FusionMatrix[gene1 + '\t' + gene2] = rawscore
            FusionPos[gene1 + '\t' + gene2] = [[chromo1, chromo2], []]
            for k in range(len(subgrp)):
                FusionPos[gene1 + '\t' + gene2][1].append([subgrp[k], AvePos[k], [i], 1])
            gene = gene1 + '\t' + gene2
    infile.close()
    exprfile.close()
'''
for key in EmcmpscoreDic:
    if key in FusionMatrix:
        for i in range(len(FusionMatrix[key])):
            FusionMatrix[key][i] += EmcmpscoreDic[key]
'''
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
        if FusionPos[gene][1][i][4][0] > -0.5:
            genepart = gene.split('\t')
            templine = genepart[0] + '\t' + genepart[1] + '\t' + str(
                FusionMatrix[gene][i]) + '\t' + FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1] + '\t' \
                       + str(FusionPos[gene][1][i][4][0]) + '\t' + str(FusionPos[gene][1][i][4][1]) + '\t' + str(i) \
                       + '\t' + str(FusionPos[gene][1][i][3])
            if EnableClustering:
                CluRes = ''
                EffClu = []
                cellset = FusionPos[gene][1][i][2]
                for cell in cellset:
                    ClusterCount[Clusters[cell]] += 1
                cclist = []
                for key in ClusterCount:
                    ClusterCount[key] /= ClusterSize[key]
                    cclist.append(ClusterCount[key])
                if max(cclist) < 0.5:
                    CluRes = "No_Spec_Clu"
                else:
                    for key in ClusterCount:
                        if ClusterCount[key] > 0.5 * max(cclist):
                            EffClu.append(key)
                    if len(EffClu) == 1:
                        CluRes = EffClu[0] + "_Spec"
                    else:
                        CluRes = "Share"
                        for item in EffClu:
                            CluRes += '_' + item
                templine += '\t' + CluRes + '\n'
            else:
                templine += '\n'
            templines.append(templine)
linedic = {}
for line in templines:
    info = line.split('\t')
    score = float(info[2])
    linedic[line] = score
for key in sorted(linedic, key=linedic.__getitem__, reverse=True):
    print(key, end='')
