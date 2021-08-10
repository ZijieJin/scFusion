from __future__ import print_function
from __future__ import division
import sys

# ##### README #######
# This code rm the genes in the fusion candidates that appears more than 3 times in the FusionScore file
# The input is the FusionScore.txt


def cmbgepos(gene, pos):
    return str(gene) + '\t' + str(pos)


FusionScoreFile = open(sys.argv[1])
genecount = {}
linestore = {}
filteredstore = {}
pospair = []
FileText = FusionScoreFile.readlines()
totalline = len(FileText)
count = 0
for line in FileText:
    if line[0] == '#':
        continue
    count += 1
    info = line.split('\t')
    gene1 = info[0]
    gene2 = info[1]
    pos1 = info[5]
    pos2 = info[6]
    rank = count / totalline
    if gene1.find(';') == -1 and gene2.find(';') == -1:
        if rank < 0.01:
            if cmbgepos(gene1, pos1) in genecount:
                genecount[cmbgepos(gene1, pos1)] += 1
            else:
                genecount[cmbgepos(gene1, pos1)] = 1
            if cmbgepos(gene2, pos2) in genecount:
                genecount[cmbgepos(gene2, pos2)] += 1
            else:
                genecount[cmbgepos(gene2, pos2)] = 1
            if cmbgepos(gene1, gene2) in pospair or cmbgepos(gene2, gene1) in pospair:
                genecount[cmbgepos(gene1, pos1)] -= 1
                genecount[cmbgepos(gene2, pos2)] -= 1
            else:
                pospair.append(cmbgepos(gene1, gene2))
        else:
            if cmbgepos(gene1, pos1) not in genecount:
                genecount[cmbgepos(gene1, pos1)] = 0
            if cmbgepos(gene2, pos2) not in genecount:
                genecount[cmbgepos(gene2, pos2)] = 0
        linestore[line] = [cmbgepos(gene1, pos1), cmbgepos(gene2, pos2)]
FusionScoreFile.close()
for key in linestore:
    if genecount[linestore[key][0]] <= 15 and genecount[linestore[key][1]] <= 15:
        filteredstore[key] = float(key.split('\t')[2])
for key in sorted(filteredstore, key=filteredstore.__getitem__, reverse=True):
    print(key, end='')
