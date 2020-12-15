from __future__ import print_function
import sys
import codecs


# ***** readme *****
# This code remove the low mappability records of Chimeric sam by star
# if the mappability < 1, delete the record locating at this position

ChimericSamFile = codecs.open(sys.argv[1], 'r',  encoding='utf-8', errors='ignore')
OutputFile = open(sys.argv[2], 'w')
MappabilityFile = open(sys.argv[3])
thres = float(sys.argv[4])
Mappability = {}
chrflag = []
keys = {}
for line in MappabilityFile.readlines():
    if line[0] == '#' or len(line) < 3:
        continue
    info = line.split('\t')
    if float(info[3].replace('\n', '').replace('\r', '')) < thres:
        continue
    if info[0] not in chrflag:
        chrflag.append(info[0])
        Mappability[info[0]] = {}
    Mappability[info[0]][int(info[1])] = int(info[2])
for i in Mappability:
    keys[i] = list(Mappability[i].keys())
    keys[i].sort()
lastname = ''
linestore = []
flag = 0
aaaaa = 0
linecount = 0
lines = ChimericSamFile.readlines()
for line in lines:
    info = line.split('\t')
    if len(info) < 5:
        continue
    if len(line) > 10:
        if line[0] == '@':
            continue
        if info[2].find('M') > -1:
            continue
    if info[0] != lastname:
        bad = 0
        poss = []
        for item in linestore:
            info = item.split('\t')
            '''
            if info[5].find('N') > -1:
                bad = 1
                break
            '''
            poss.append(info[3])
        if len(poss) <= 1 or (len(poss) == 2 and abs(int(poss[0])-int(poss[1])) <= 10):
            bad = 1
        if bad == 0:
            mapflag = 1
            index = 0
            for subline in linestore:
                if (index == 1 and poss[0] == poss[1]) or (index == 2 and (poss[2] == poss[1] or poss[2] == poss[0])):
                    index += 1
                    continue
                index += 1
                info = subline.split('\t')
                if not info[2].startswith('chr'):
                    info[2] = 'chr' + info[2]
                findstart = 0
                if info[2] in keys:
                    findend = len(keys[info[2]]) - 1
                    if keys[info[2]][findstart] > int(info[3]):
                        mapflag = 0
                        break
                    if Mappability[info[2]][keys[info[2]][findend]] <= int(info[3]):
                        mapflag = 0
                        break
                    if keys[info[2]][findstart] <= int(info[3]) <= Mappability[info[2]][keys[info[2]][findstart]]:
                        continue
                    if keys[info[2]][findend] <= int(info[3]) <= Mappability[info[2]][keys[info[2]][findend]]:
                        continue
                    whilecount = 0
                    while True:
                        whilecount += 1
                        if findend - findstart <= 1:
                            mapflag = 0
                            break
                        mid = int((findend - findstart) * 0.618) + findstart
                        if keys[info[2]][mid] > int(info[3]):
                            findend = mid
                        elif keys[info[2]][mid] <= int(info[3]):
                            if Mappability[info[2]][keys[info[2]][mid]] <= int(info[3]):
                                findstart = mid
                            else:
                                ppp = mid
                                break
            if mapflag == 1:
                for subline in linestore:
                    flag += 1
                    try:
                        sep = '\t'
                        subinfo = subline.split('\t')
                        if not subinfo[2].startswith('chr'):
                            subinfo[2] = 'chr' + subinfo[2]
                        OutputFile.write(sep.join(subinfo))
                    except:
                        print(subline)
        info = line.split('\t')
        lastname = info[0]
        linestore = [line]
    else:
        linestore.append(line)
    if line == '':
        break
ChimericSamFile.close()
MappabilityFile.close()
OutputFile.close()
