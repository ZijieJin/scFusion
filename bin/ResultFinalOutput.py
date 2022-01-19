from __future__ import print_function
from __future__ import division
import sys

# ***** readme *****
# 添加信息： read，direct，insertsize, anchor

infile = open(sys.argv[1])
refannotfile = open(sys.argv[2])
ChimericOutDir = sys.argv[3]
OutPrefix = sys.argv[4]
#nodetailcell = sys.argv[5]
nodetailcell = 0
outfull = open(OutPrefix + '.full.txt', 'w')
outabridge = open(OutPrefix + '.abridged.txt', 'w')
genestrand = {}
genecount = {}
repeatcount = {}
posrecord = []
for line in refannotfile.readlines():
    if line.startswith('#'):
        continue
    info = line.split('\t')
    if info[2] == 'exon':
        chromo = info[0]
        start = int(info[3])
        end = int(info[4])
        strand = info[6]
        geneinfo = info[8].split('; ')
        genename = ''
        for item in geneinfo:
            if item.startswith('gene_name'):
                genename = item[11:-1]
            if item.startswith('gene_type'):
                genetype = item[11:-1]
            if item.startswith('gene_biotype'):
                genetype = item[14:-1]
        if genename == '':
            for item in geneinfo:
                if item.find('gene_id') > -1:
                    genename = item[9:-1]
        if genename == '':
            continue
        genestrand[genename] = strand
refannotfile.close()
outabridge.write(
    '#Fusion_id\tFusionGene\tPosition1\tPosition2\tstrands\tdirections\tSupportingCellNumber\tTotalSplitRead\tTotalDiscordantRead\tFakeProbability\tFDR\n')
outfull.write(
    '#Fusion_id\tFusionGene\tPosition1\tPosition2\tstrands\tdirections\tSupportingCellNumber\tSupportingCells\tTotalSplitRead\tSplitReads\tTotalDiscordantRead\tDiscordantReads\tFakeProbability\tFDR\n')
currentid = 0
lines = infile.readlines()
uselines = []
finaluselines = []
for line in lines:
    if line.startswith('#'):
        continue
    resinfo = line.rstrip('\n').split('\t')
    fusiongene1 = resinfo[0].split('--')[0]
    fusiongene2 = resinfo[0].split('--')[1]
    fusionpos1 = int(resinfo[4].split(':')[1])
    fusionpos2 = int(resinfo[5].split(':')[1])
    '''
    found = False
    for item in posrecord:
        if abs(fusionpos1 - item[0]) == abs(fusionpos2 - item[1]) or abs(fusionpos2 - item[0]) == abs(
                fusionpos1 - item[1]) or fusionpos2 + item[1] == fusionpos1 + item[0] or fusionpos1 + item[
            1] == fusionpos2 + item[0]:
            found = True
            break
    if not found:
        posrecord.append([fusionpos1, fusionpos2])
        uselines.append(line)
        if fusiongene1 not in genecount:
            genecount[fusiongene1] = 0
        genecount[fusiongene1] += 1
        if fusiongene2 not in genecount:
            genecount[fusiongene2] = 0
        genecount[fusiongene2] += 1
    '''
    if fusiongene1 not in repeatcount:
        repeatcount[fusiongene1] = {}
    if fusiongene2 not in repeatcount:
        repeatcount[fusiongene2] = {}
    if fusionpos1 not in repeatcount[fusiongene1]:
        repeatcount[fusiongene1][fusionpos1] = []
    if fusionpos2 not in repeatcount[fusiongene2]:
        repeatcount[fusiongene2][fusionpos2] = []
    if fusiongene2 not in repeatcount[fusiongene1][fusionpos1]:
        repeatcount[fusiongene1][fusionpos1].append(fusiongene2)
    if fusiongene1 not in repeatcount[fusiongene2][fusionpos2]:
        repeatcount[fusiongene2][fusionpos2].append(fusiongene1)
'''
for line in uselines:
    if line.startswith('#'):
        continue
    resinfo = line.rstrip('\n').split('\t')
    fusiongene1 = resinfo[0].split('--')[0]
    fusiongene2 = resinfo[0].split('--')[1]
    if genecount[fusiongene1] >= 5 or genecount[fusiongene2] >= 5:
        continue
    finaluselines.append(line)
'''
for line in lines:
    if line.startswith('#'):
        continue
    resinfo = line.rstrip('\n').split('\t')
    cellsup = resinfo[-1].split(', ')[:-1]
    fusiongene1 = resinfo[0].split('--')[0]
    fusiongene2 = resinfo[0].split('--')[1]
    fusionpos1 = int(resinfo[4].split(':')[1])
    fusionpos2 = int(resinfo[5].split(':')[1])
    direct1 = resinfo[9]
    direct2 = resinfo[10]
    try:
        strand1 = genestrand[fusiongene1]
        if direct1 == strand1:
            direct1 = 'd'
        else:
            direct1 = 'u'
    except KeyError:
        strand1 = 'N/A'
        direct1 = 'N/A'
    try:
        strand2 = genestrand[fusiongene2]
        if direct2 == strand2:
            direct2 = 'd'
        else:
            direct2 = 'u'
    except KeyError:
        strand2 = 'N/A'
        direct2 = 'N/A'
    anchor = []
    insertsize = []
    discordant = []
    splitreadnum = []
    totalsplitread = 0
    totaldiscordant = 0
    if len(repeatcount[fusiongene1][fusionpos1]) >= 5 or len(repeatcount[fusiongene2][fusionpos2]) >= 5:
        continue
    if nodetailcell != '1':
        for cellsupitem in cellsup:
            if cellsupitem == '' or cellsupitem == ' ':
                continue
            try:
                fsfile = open(ChimericOutDir + '/' + cellsupitem + '_FusionSupport.txt')
                for fsline in fsfile.readlines():
                    info = fsline.rstrip().split('\t')
                    gene1 = info[0]
                    gene2 = info[1]
                    if gene1 == fusiongene1 and gene2 == fusiongene2 or gene1 == fusiongene2 and gene2 == fusiongene1:
                        discordant.append(int(info[2]))
                        if len(info) < 7:
                            splitsup = []
                        else:
                            splitsup = info[6].split(';')
                        splitcount = 0
                        for item in splitsup:
                            if item.find('+') == -1:
                                continue
                            thispos1 = int(item.split('+')[0].split(',')[0])
                            thispos2 = int(item.split('+')[0].split(',')[1])
                            if abs(thispos2 - fusionpos2) + abs(thispos1 - fusionpos1) < 30 or abs(
                                    thispos2 - fusionpos1) + abs(thispos1 - fusionpos2) < 30:
                                splitcount += 1
                        splitreadnum.append(splitcount)
                        break
            except IOError:
                pass
        totalsplitread = sum(splitreadnum)
        totaldiscordant = sum(discordant)
        cellsup = list(map(int, cellsup))
    currentid += 1
    if direct1 == 'd' and direct2 == 'u':
        fgenes = resinfo[0].split('--')
        outabridge.write(str(currentid) + '\t' + fgenes[1] + '--' + fgenes[0] + '\t' + resinfo[5] + '\t' + resinfo[
            4] + '\t' + strand2 + '/' + strand1 + '\t' + direct2 + '/' + direct1 + '\t' +
                         str(len(splitreadnum)) + '\t' + str(totalsplitread) + '\t' + str(totaldiscordant) + '\t' + resinfo[
                             6] + '\t' + resinfo[8] + '\n')
        outfull.write(str(currentid) + '\t' + fgenes[1] + '--' + fgenes[0] + '\t' + resinfo[5] + '\t' + resinfo[
            4] + '\t' + strand2 + '/' + strand1 + '\t' + direct2 + '/' + direct1 + '\t' + str(
            len(splitreadnum)) + '\t' + str(cellsup) + '\t' + str(totalsplitread) + '\t' + str(splitreadnum) + '\t' +
                      str(totaldiscordant) + '\t' + str(discordant) + '\t' + resinfo[6] + '\t' + resinfo[8] + '\n')
    else:
        outabridge.write(str(currentid) + '\t' + resinfo[0] + '\t' + resinfo[4] + '\t' + resinfo[
            5] + '\t' + strand1 + '/' + strand2 + '\t' + direct1 + '/' + direct2 + '\t' +
                         str(len(splitreadnum)) + '\t' + str(totalsplitread) + '\t' + str(totaldiscordant) + '\t' + resinfo[
                             6] + '\t' + resinfo[8] + '\n')
        outfull.write(str(currentid) + '\t' + resinfo[0] + '\t' + resinfo[4] + '\t' + resinfo[
            5] + '\t' + strand1 + '/' + strand2 + '\t' + direct1 + '/' + direct2 + '\t' + str(
            len(splitreadnum)) + '\t' + str(cellsup) + '\t' + str(totalsplitread) + '\t' + str(splitreadnum) + '\t' +
                      str(totaldiscordant) + '\t' + str(discordant) + '\t' + resinfo[6] + '\t' + resinfo[8] + '\n')
infile.close()
outabridge.close()
outfull.close()
