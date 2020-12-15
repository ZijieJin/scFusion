from __future__ import print_function
import sys
# ***** readme *****
# This code analyses softclipping to find encompassing reads and split reads for one sam file




def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


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


def TakeoutFusionSupport(lines):
    if len(lines) == 2:
        info1 = lines[0].split('\t')
        info2 = lines[1].split('\t')
        if info1[0] != info2[0]:
            cgene1 = info1[0]
            cgene2 = info2[0]
            chromo1 = info1[3]
            chromo2 = info2[3]
            if not chromo1.startswith('chr'):
                chromo1 = 'chr' + chromo1
            if not chromo2.startswith('chr'):
                chromo2 = 'chr' + chromo2
            genelist1 = cgene1.split(';')
            genelist2 = cgene2.split(';')  # separate genes with semicolon
            for gene1 in genelist1:
                for gene2 in genelist2:
                    if gene1 + '\t' + gene2 in geneset:
                        geneset[gene1 + '\t' + gene2][0] += 1
                    elif gene2 + '\t' + gene1 in geneset:
                        geneset[gene2 + '\t' + gene1][0] += 1
                    else:
                        geneset[gene1 + '\t' + gene2] = [1, 0, [], [chromo1, chromo2]]
    elif len(lines) == 3:
        info1 = lines[0].split('\t')
        info2 = lines[1].split('\t')
        info3 = lines[2].split('\t')
        gene1 = info1[0]
        gene2 = info2[0]
        gene3 = info3[0]
        if (gene1 != gene2 or gene2 != gene3) and (gene1 == gene2 or gene2 == gene3 or gene1 == gene3) and (gene1 != '' and gene2 != '' and gene3 != ''):
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
                print('Three reads are different!')
                for ll in lines:
                    print(ll, end='')
                return thisname
            if readinfo1[0] != readinfo2[0]:
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
                if readlen - 0.5 > splitpnt1 > -0.5 and readlen - 0.5 > splitpnt2 > 0.5:
                    if clipsplit1[1][i] == 'S':
                        brkpnt1 = int(readinfo1[4])
                        direct1 = 1
                    else:
                        brkpnt1 = splitpnt1 - clipsplit1[2] + int(readinfo1[4]) + clipsplit1[4] - 1
                        direct1 = -1
                    if clipsplit2[1][j] == 'S':
                        brkpnt2 = int(readinfo2[4])
                        direct2 = 1
                    else:
                        brkpnt2 = splitpnt2 - clipsplit2[2] + int(readinfo2[4]) + clipsplit2[4] - 1
                        direct2 = -1
                    cgene1 = readinfo1[0]
                    cgene2 = readinfo2[0]
                    chromo1 = readinfo1[3]
                    chromo2 = readinfo2[3]
                    if not chromo1.startswith('chr'):
                        chromo1 = 'chr' + chromo1
                    if not chromo2.startswith('chr'):
                        chromo2 = 'chr' + chromo2
                    genelist1 = cgene1.split(';')
                    genelist2 = cgene2.split(';')  # separate genes with semicolon
                    for gene1 in genelist1:
                        for gene2 in genelist2:
                            if gene1 + '\t' + gene2 in geneset:
                                geneset[gene1 + '\t' + gene2][1] += 1
                                geneset[gene1 + '\t' + gene2][2].append(
                                    [brkpnt1, brkpnt2, splitpnt1, clipsplit1[2], clipsplit2[2], direct1, direct2])
                            elif gene2 + '\t' + gene1 in geneset:
                                geneset[gene2 + '\t' + gene1][1] += 1
                                geneset[gene2 + '\t' + gene1][2].append(
                                    [brkpnt2, brkpnt1, splitpnt2, clipsplit2[2], clipsplit1[2], direct2, direct1])
                            else:
                                geneset[gene1 + '\t' + gene2] = [0, 1,
                                                                 [[brkpnt1, brkpnt2, splitpnt1, clipsplit1[2], clipsplit2[2], direct1, direct2]],
                                                                 [chromo1, chromo2]]
    return ''


geneset = {}
infile = open(sys.argv[1])
outfile = open(sys.argv[2], 'w')
lastname = ''
lines = []
badreadname = ''
thisname = 'hgfhfrjfjzjbest1122gh'
with infile:
    for line in infile:
        if len(line) > 0:
            if line[0] == '#':
                continue
        info = line.split('\t')
        if len(line) > 0:
            thisname = info[1]
        if thisname == badreadname:
            continue
        if thisname == lastname:
            lines.append(line)
        else:
            aa = TakeoutFusionSupport(lines)  # mean procedure to record fusion.
            if len(aa) > 1:
                lastname = aa
            else:
                lastname = thisname
            lines = [line]
        if len(info[0]) <= 1:
            badreadname = info[1]
            lastname = ''
            lines = []
aa = TakeoutFusionSupport(lines)
for key in geneset:
    outfile.write(key + '\t' + str(geneset[key][0]) + '\t' + str(geneset[key][1]) + '\t' + geneset[key][3][0] + '\t' + geneset[key][3][1] + '\t')
    posstr = []
    for i in geneset[key][2]:
        # posstr.append(str(i[0]) + ',' + str(i[1]) + '+' + str(i[2]) + '+' + str(i[3]) + ',' + str(i[4]) + '+' + str(i[5]) + ',' + str(i[6]) + ';')
        posstr.append(str(i[0]) + ',' + str(i[1]) + '+' + str(i[2]) + '+' + str(i[3]) + ',' + str(i[4]) + ';')
    posstr = set(posstr)
    for i in posstr:
        outfile.write(i)
    outfile.write('\n')
infile.close()
outfile.close()
