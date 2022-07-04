from __future__ import print_function
import sys
import codecs


# ***** readme *****
# This code does not remove the low mappability records of Chimeric sam by star

ChimericSamFile = codecs.open(sys.argv[1], 'r',  encoding='utf-8', errors='ignore')
OutputFile = open(sys.argv[2], 'w')
chrflag = []
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
    if not info[2].startswith('chr'):
        info[2] = 'chr' + info[2]
    OutputFile.write(line)
    if line == '':
        break
ChimericSamFile.close()
OutputFile.close()
