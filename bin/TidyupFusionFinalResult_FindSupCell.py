from __future__ import print_function
from __future__ import division
import sys
import os.path


# ***** readme *****
# This code removes the duplicates and add supporting info.

def CheckGoodGene(gene, pos):
    if not FilteringSwitch:
        return True
    sum = 0
    for item in FusionGeneRecord[gene]:
        sum += item[1]
    for item in FusionGeneRecord[gene]:
        if abs(item[0] - pos) < 5:
            if item[1] / sum < 0.4 and item[2] > 3:
                return False
            return True
    return False

FilteringSwitch = True
ResultFile = open(sys.argv[1])
if len(sys.argv) > 3:
    FilteringSwitch = False
FileDir = sys.argv[2]
PosRecord = []
FusionGeneRecord = {}
for line in ResultFile.readlines():
    if line[0] == '#':
        continue
    if line.startswith('FusionName'):
        continue
    info = line.split('\t')
    pos1 = int(info[4].split(':')[1])
    pos2 = int(info[5].split(':')[1].rstrip('\n'))
    gene = info[0].split('--')
    gene1 = gene[0]
    gene2 = gene[1]
    if [info[5], info[4]] in PosRecord or [info[4], info[5]] in PosRecord:
        continue
    PosRecord.append([info[4], info[5]])
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
ResultFile = open(sys.argv[1])
PosRecord = []
for line in ResultFile.readlines():
    if line[0] == '#':
        print(line[:-1] + '\tSupportingCells')
        continue
    if line.startswith('FusionName'):
        continue
    info = line.split('\t')
    if [info[5], info[4]] in PosRecord or [info[4], info[5]] in PosRecord:
        continue
    PosRecord.append([info[4], info[5]])
    gene = info[0].split('--')
    gene1 = gene[0]
    gene2 = gene[1]
    pos1 = int(info[4].split(':')[1])
    pos2 = int(info[5].split(':')[1].rstrip('\n'))
    if CheckGoodGene(gene1, pos1) and CheckGoodGene(gene2, pos2):
        print(line[:-1], end='\t')
        cellname = []
        filenames = os.listdir(FileDir)
        for file in filenames:
            if file.endswith('_FusionSupport.txt'):
                thisfile = open(FileDir + file)
                for thisline in thisfile.readlines():
                    info = thisline.split('\t')
                    found = False
                    if len(info) == 7 and len(info[6]) > 3 and ((gene1 == info[0] and gene2 == info[1]) or (gene2 == info[0] and gene1 == info[1])):
                        readsup = info[6].rstrip(';')
                        splitinfo = readsup.split(';')
                        for item in splitinfo:
                            if len(item) <= 3:
                                continue
                            a = item.split('+')[0]
                            pp = a.split(',')
                            try:
                                if abs(pos1 - int(pp[0])) < 10 and abs(pos2 - int(pp[1])) < 10 or abs(pos2 - int(pp[0])) < 10 and abs(pos1 - int(pp[1])) < 10:
                                    cellname.append(int(file.split('_')[0]))
                                    found = True
                                    break
                            except:
                                sys.stderr.write(str(pos1) + '\t' + str(pos2) + '\t' + item)
                    elif len(info) == 6:
                        if gene1 == info[0] and gene2 == info[1] or gene2 == info[0] and gene1 == info[1]:
                            cellname.append(int(file.split('_')[0]))
                    if found:
                        break
                thisfile.close()
        cellname = sorted(cellname)
        for item in cellname:
            print(str(item), end=', ')
        print('')
ResultFile.close()
