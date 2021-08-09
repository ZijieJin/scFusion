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
FindSupCell = True
ResultFile = open(sys.argv[1])
if len(sys.argv) > 3:
    FilteringSwitch = False
FileDir = sys.argv[2]
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


FusionRecCells = {}
filecount = 0
filenames = os.listdir(FileDir)
for file in filenames:
    if file.endswith('_FusionSupport.txt'):
        thisfile = open(FileDir + file)
        filecount += 1
        cellindex = int(file.split('_')[0])
        FusionRecCells[cellindex] = {}
        for thisline in thisfile.readlines():
            info = thisline.split('\t')
            found = False
            if len(info) == 7:
                readsup = info[6].rstrip(';')
                splitinfo = readsup.split(';')
                FusionRecCells[cellindex][info[0] + '--' + info[1]] = splitinfo
            elif len(info) == 6:
                FusionRecCells[cellindex][info[0] + '--' + info[1]] = []
        thisfile.close()



for line in uselines:
    if line[0] == '#':
        continue
    if line.startswith('FusionName'):
        continue
    info = line.split('\t')
    gene = info[0].split('--')
    gene1 = gene[0]
    gene2 = gene[1]
    pos1 = int(info[4].split(':')[1])
    pos2 = int(info[5].split(':')[1].rstrip('\n'))
    if CheckGoodGene(gene1, pos1) and CheckGoodGene(gene2, pos2):
        print(line[:-1], end='\t')
        cellname = []
        for cell in FusionRecCells:
            if gene1 + '--' + gene2 in FusionRecCells[cell]:
                key = gene1 + '--' + gene2
            elif gene2 + '--' + gene1 in FusionRecCells[cell]:
                key = gene2 + '--' + gene1
            else:
                continue
            if FusionRecCells[cell][key] == []:
                cellname.append(cell)
            else:
                for item in FusionRecCells[cell][key]:
                    if len(item) <= 3:
                        continue
                    a = item.split('+')[0]
                    pp = a.split(',')
                    try:
                        if abs(pos1 - int(pp[0])) < 10 and abs(pos2 - int(pp[1])) < 10 or abs(pos2 - int(pp[0])) < 10 and abs(pos1 - int(pp[1])) < 10:
                            cellname.append(cell)
                            break
                    except:
                        sys.stderr.write(str(pos1) + '\t' + str(pos2) + '\t' + item)
        cellname = sorted(cellname)
        if cellname == []:
            print('NoSupCell', end='')
        else:
            for item in cellname:
                print(str(item), end=', ')
        print('')
ResultFile.close()
