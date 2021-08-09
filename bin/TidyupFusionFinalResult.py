from __future__ import print_function
from __future__ import division
import sys
import os.path


# ***** readme *****
# This code removes the duplicates and add supporting info.

def CheckGoodGene(gene, pos, thiscount):
    if not FilteringSwitch:
        return True
    sum = 0
    for item in FusionGeneRecord[gene]:
        sum += item[1]
    for item in FusionGeneRecord[gene]:
        if abs(item[0] - pos) < 5:
            if thiscount / sum < 0.4 and item[2] >= 3:
                return False
            return True
    return False


FilteringSwitch = True
ResultFile = open(sys.argv[1])
if len(sys.argv) > 2:
    FilteringSwitch = False
PosRecord = {}
FusionGeneRecord = {}
uselines = []
for line in ResultFile.readlines():
    if line[0] == '#':
        print(line[:-1] + '\tSupportingCells')
        continue
    if line.startswith('FusionName'):
        print(line[:-1] + '\tSupportingCells')
        continue
    info = line.split('\t')
    pos1 = int(info[4].split(':')[1])
    pos2 = int(info[5].split(':')[1].rstrip('\n'))
    gene = info[0].split('--')
    gene1 = gene[0]
    gene2 = gene[1]
    '''
    if info[5] + '--' + info[4] in PosRecord or info[4] + '--' + info[5] in PosRecord:
        continue
    '''
    PosRecord[info[4] + '--' + info[5]] = 0
    uselines.append(line)
    if gene1 not in FusionGeneRecord:
        FusionGeneRecord[gene1] = [[pos1, int(info[2]), 1]]
    else:
        near = False
        for i in range(len(FusionGeneRecord[gene1])):
            if abs(FusionGeneRecord[gene1][i][0] - pos1) <= 5:
                near = True
                FusionGeneRecord[gene1][i][1] += int(info[2])
                FusionGeneRecord[gene1][i][2] += 1
                break
        if not near:
            FusionGeneRecord[gene1].append([pos1, int(info[2]), 1])
    if gene2 not in FusionGeneRecord:
        FusionGeneRecord[gene2] = [[pos2, int(info[2]), 1]]
    else:
        near = False
        for i in range(len(FusionGeneRecord[gene2])):
            if abs(FusionGeneRecord[gene2][i][0] - pos2) <= 5:
                near = True
                FusionGeneRecord[gene2][i][1] += int(info[2])
                FusionGeneRecord[gene2][i][2] += 1
                break
        if not near:
            FusionGeneRecord[gene2].append([pos2, int(info[2]), 1])
ResultFile.close()
PosRecord = []
for line in uselines:
    if line[0] == '#':
        continue
    if line.startswith('FusionName'):
        continue
    info = line.split('\t')
    gene = info[0].split('--')
    gene1 = gene[0]
    gene2 = gene[1]
    gene1main = gene1.split('-')[0]
    gene2main = gene2.split('-')[0]
    pos1 = int(info[4].split(':')[1])
    pos2 = int(info[5].split(':')[1].rstrip('\n'))
    if CheckGoodGene(gene1, pos1, int(info[2])) and CheckGoodGene(gene2, pos2, int(info[2])) and gene1main != gene2main:
        print(line[:-1])
ResultFile.close()
