from __future__ import print_function
from __future__ import division
import sys
import math


# ***** readme *****
# This code add transcriptome read to the simulated file


def findpos(chromo, pos):
    row = int(math.ceil(pos / reffilelinelength))
    col = pos % reffilelinelength
    if col == 0:
        col = reffilelinelength
    return row, col


def getTCGAstring(chr, start, end):
    row1, col1 = findpos(chr, start)
    row2, col2 = findpos(chr, end)
    if row1 == row2:
        return refdict[chr][row1 - 1][col1 - 1:col2]
    resstring = refdict[chr][row1 - 1].rstrip()[col1 - 1:]
    for i in range(row1 + 1, row2):
        resstring += refdict[chr][i - 1].rstrip()
    return resstring + refdict[chr][row2 - 1][:col2]


rownum = [1, 4154180, 8207504, 11507879, 14693785, 17709041, 20560960, 23213273, 25652675, 28006234, 30265148, 32515258,
          34746124, 36665623, 38454783, 40163641, 41669555, 43022810, 44324099, 45309584, 46360011, 47162177, 48017255,
          50605099]

refannotfile = open(sys.argv[1])
hg19file = open(sys.argv[2])
readsfile = open(sys.argv[3])
halfwindowslength = 300
geneset = {}
PosRecord = {}

for line in refannotfile.readlines():
    if line.startswith('#'):
        continue
    info = line.split('\t')
    if info[2] == 'exon':
        chromo = info[0]
        start = int(info[3])
        end = int(info[4])
        geneinfo = info[8].split('; ')
        genename = ''
        if len(geneinfo) > 4:
            for item in geneinfo:
                if item.find('gene_name') > -1:
                    genename = item[11:-1]
                if item.find('gene_type') > -1:
                    genetype = item[11:-1]
                if item.find('gene_biotype') > -1:
                    genetype = item[14:-1]
            if genename == '':
                for item in geneinfo:
                    if item.find('gene_id') > -1:
                        genename = item[9:-1]
        else:
            continue
        if genetype != 'protein_coding':
            continue
        if chromo not in geneset:
            geneset[chromo] = {}
        if genename == '':
            continue
        if genename not in geneset[chromo]:
            geneset[chromo][genename] = [[start, end], [[start, end]]]
        else:
            geneset[chromo][genename][1].append([start, end])
            geneset[chromo][genename][0][0] = min(geneset[chromo][genename][0][0], start)
            geneset[chromo][genename][0][1] = max(geneset[chromo][genename][0][1], end)
refannotfile.close()

reference = hg19file.readlines()

currentchr = ''
refdict = {}
for line in reference:
    if line.startswith('>'):
        currentchr = line.rstrip()[1:].split(' ')[0]
        if not currentchr.startswith('chr'):
            currentchr = 'chr' + currentchr
        refdict[currentchr] = []
    else:
        refdict[currentchr].append(line)
reffilelinelength = len(refdict[currentchr][0].rstrip())
if reffilelinelength > 100:
    reffilelinelength = 99999999999


for chromo in geneset:
    for genename in geneset[chromo]:
        for i in range(len(geneset[chromo][genename][1])):
            for j in range(i+1, len(geneset[chromo][genename][1])):
                if geneset[chromo][genename][1][j][0] < geneset[chromo][genename][1][i][0]:
                    geneset[chromo][genename][1][j][0], geneset[chromo][genename][1][i][0] = geneset[chromo][genename][1][i][0], geneset[chromo][genename][1][j][0]
                    geneset[chromo][genename][1][j][1], geneset[chromo][genename][1][i][1] = geneset[chromo][genename][1][i][1], geneset[chromo][genename][1][j][1]
        k = 0
        while k < len(geneset[chromo][genename][1]) - 1:
            if geneset[chromo][genename][1][k + 1][0] <= geneset[chromo][genename][1][k][1]:
                geneset[chromo][genename][1][k][1] = geneset[chromo][genename][1][k + 1][1]
                del geneset[chromo][genename][1][k+1]
            else:
                k += 1



for line in readsfile.readlines():
    info = line.split('\t')
    chr1 = info[2].split(":")[0]
    pos1 = int(info[2].split(':')[1])
    chr2 = info[3].split(":")[0]
    pos2 = int(info[3].split(':')[1])
    gene1 = ''
    gene2 = ''
    exoninterval1 = -1
    exoninterval2 = -1
    if chr1 not in geneset or chr2 not in geneset:
        continue
    for g1 in geneset[chr1]:
        if geneset[chr1][g1][0][0] <= pos1 <= geneset[chr1][g1][0][1]:
            gene1 = g1
            break
    for g2 in geneset[chr2]:
        if geneset[chr2][g2][0][0] <= pos2 <= geneset[chr2][g2][0][1]:
            gene2 = g2
            break
    if gene1 != '' and gene2 != '':
        for i in range(len(geneset[chr1][gene1][1])):
            if geneset[chr1][gene1][1][i][0] <= pos1 <= geneset[chr1][gene1][1][i][1]:
                exoninterval1 = i
                break
        for i in range(len(geneset[chr2][gene2][1])):
            if geneset[chr2][gene2][1][i][0] <= pos2 <= geneset[chr2][gene2][1][i][1]:
                exoninterval2 = i
                break
    if exoninterval2 == -1 or exoninterval1 == -1:
        res1d = getTCGAstring(chr1, pos1, pos1 + halfwindowslength)
        res1u = getTCGAstring(chr1, pos1 - halfwindowslength, pos1-1)
        res2d = getTCGAstring(chr2, pos2, pos2 + halfwindowslength)
        res2u = getTCGAstring(chr2, pos2 - halfwindowslength, pos2-1)
    else:
        res1d = ''
        res2d = ''
        res1u = ''
        res2u = ''
        last1 = 0
        last2 = 0
        if geneset[chr1][gene1][1][exoninterval1][1] - pos1 + 1 >= halfwindowslength + 1:
            res1d = getTCGAstring(chr1, pos1, pos1+halfwindowslength-1)
        else:
            res1d = getTCGAstring(chr1, pos1, geneset[chr1][gene1][1][exoninterval1][1])
            last1 = geneset[chr1][gene1][1][exoninterval1][1]
        if geneset[chr2][gene2][1][exoninterval2][1] - pos1 + 1 >= halfwindowslength + 1:
            res2d = getTCGAstring(chr2, pos2, pos2+halfwindowslength-1)
        else:
            res2d = getTCGAstring(chr2, pos2, geneset[chr2][gene2][1][exoninterval2][1])
            last2 = geneset[chr2][gene2][1][exoninterval2][1]
        fi1 = exoninterval1
        fi2 = exoninterval2
        while len(res1d) < halfwindowslength:
            fi1 += 1
            if fi1 >= len(geneset[chr1][gene1][1]):
                break
            if last1 >= geneset[chr1][gene1][1][fi1][0]:
                break
            if geneset[chr1][gene1][1][fi1][1] - geneset[chr1][gene1][1][fi1][0] + 1 + len(res1d) >= halfwindowslength + 1:
                res1d += getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0], geneset[chr1][gene1][1][fi1][0] + halfwindowslength - len(res1d) - 1)
            else:
                res1d += getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0], geneset[chr1][gene1][1][fi1][1])
                last1 = geneset[chr1][gene1][1][fi1][1]

        if len(res1d) < halfwindowslength:
            res1d += getTCGAstring(chr1, last1 + 1, last1 + halfwindowslength - len(res1d))

        while len(res2d) < halfwindowslength:
            fi2 += 1
            if fi2 >= len(geneset[chr2][gene2][1]):
                break
            if last2 >= geneset[chr2][gene2][1][fi2][0]:
                break
            if geneset[chr2][gene2][1][fi2][1] - geneset[chr2][gene2][1][fi2][0] + 1 + len(res2d) >= halfwindowslength + 1:
                res2d += getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0], geneset[chr2][gene2][1][fi2][0] + halfwindowslength - len(res2d) - 1)
            else:
                res2d += getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0], geneset[chr2][gene2][1][fi2][1])
                last2 = geneset[chr2][gene2][1][fi2][1]

        if len(res2d) < halfwindowslength:
            res2d += getTCGAstring(chr2, last2 + 1, last2 + halfwindowslength - len(res2d))

        if - geneset[chr1][gene1][1][exoninterval1][0] + pos1 >= halfwindowslength:
            res1u = getTCGAstring(chr1, pos1-halfwindowslength, pos1-1)
        else:
            res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][exoninterval1][0], pos1 - 1)
            last1 = geneset[chr1][gene1][1][exoninterval1][0]
        if - geneset[chr2][gene2][1][exoninterval2][0] + pos2 >= halfwindowslength:
            res2u = getTCGAstring(chr2, pos2-halfwindowslength, pos2-1)
        else:
            res2u = getTCGAstring(chr2, geneset[chr2][gene2][1][exoninterval2][0], pos2 - 1)
            last2 = geneset[chr2][gene2][1][exoninterval2][0]
        fi1 = exoninterval1
        fi2 = exoninterval2
        while len(res1u) < halfwindowslength:
            fi1 -= 1
            if fi1 < 0:
                break
            if last1 <= geneset[chr1][gene1][1][fi1][1]:
                break
            if geneset[chr1][gene1][1][fi1][1] - geneset[chr1][gene1][1][fi1][0] + 1 + len(res1u) >= halfwindowslength:
                res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][1] - halfwindowslength + len(res1u) + 1, geneset[chr1][gene1][1][fi1][1]) + res1u
            else:
                res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0], geneset[chr1][gene1][1][fi1][1]) + res1u
                last1 = geneset[chr1][gene1][1][fi1][0]

        if len(res1u) < halfwindowslength:
            res1u = getTCGAstring(chr1, last1 - halfwindowslength + len(res1u), last1 - 1) + res1u

        while len(res2u) < halfwindowslength:
            fi2 -= 1
            if fi2 < 0:
                break
            if last2 <= geneset[chr2][gene2][1][fi2][1]:
                break
            if geneset[chr2][gene2][1][fi2][1] - geneset[chr2][gene2][1][fi2][0] + 1 + len(res2u) >= halfwindowslength:
                res2u = getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][1] - halfwindowslength + len(res2u) + 1, geneset[chr2][gene2][1][fi2][1]) + res2u
            else:
                res2u = getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0], geneset[chr2][gene2][1][fi2][1]) + res2u
                last2 = geneset[chr2][gene2][1][fi2][0]

        if len(res2u) < halfwindowslength:
            res2u = getTCGAstring(chr2, last2 - halfwindowslength + len(res2u), last2 - 1) + res2u

    res1 = res1u.lower() + res1d.upper()
    res2 = res2u.lower() + res2d.upper()
    print(line.rstrip() + '\t' + res1 + '\t' + res2)
readsfile.close()
