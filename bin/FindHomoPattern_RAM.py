from __future__ import print_function
from __future__ import division
import sys
import numpy as np
import math
import random

# ***** readme *****
# This code wants to verify that fake fusion is more likely to be homo


def simulateString(length):
    baseset = ["A", "T", "C", "G"]
    res = ''
    for i in range(length):
        res += str(random.sample(baseset, 1))[2]
    return res


def GCcontent(str1, str2):
    l1 = len(str1)
    l2 = len(str2)
    count = 0
    for i in range(l1):
        if str1[i].upper() == 'G' or str1[i].upper() == 'C':
            count += 1
    res1 = count / l1
    count = 0
    for i in range(l2):
        if str2[i].upper() == 'G' or str2[i].upper() == 'C':
            count += 1
    res2 = count / l2
    return [res1, res2]


def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


def Findbrkpnt_SWResult(string, tailpos, halfwidth):
    aa = []
    for i in range(len(string)):
        if string[i] != '-':
            aa.append(string[i])
    if tailpos < halfwidth:
        return len(string) + halfwidth - tailpos
    if tailpos - len(aa) >= halfwidth:
        return - tailpos + len(aa) + halfwidth
    j = len(string) - 1
    start = tailpos
    while tailpos > halfwidth:
        if string[j] != '-':
            start -= 1
        j -= 1
        if start == halfwidth:
            break
    return j + 1

'''
def myscore(str1, str2):    # added pos-by-pos
    sc = 0
    if len(str1) * len(str2) == 0:
        return 0
    for i in range(min(len(str1), len(str2))):
        if str1[i] == str2[i]:
            sc += 1
    return sc


def myscore2(str1, str2, mismatch=1):   # find a longest same substring with tolerance mismatch
    sc = 0
    if len(str1) * len(str2) == 0:
        return 0
    for i in range(min(len(str1), len(str2))):
        sctemp = 0
        miscount = 0
        for j in range(i, min(len(str1), len(str2))):
            badflag = 0
            if str1[j] == str2[j]:
                sctemp += 1
                badflag = 0
            else:
                miscount += 1
                if miscount > mismatch:
                    sctemp -= badflag * 0.1
                    break
                sctemp += 0.1
                badflag += 1
        if sc < sctemp:
            sc = sctemp
    return sc


def myscore3(str1, str2):
    len1 = len(str1)
    midpoint = math.ceil(len1 / 2) - 0.5
    score = []
    if len(str1) * len(str2) == 0:
        return 0
    for i in range(len1):
        if str1[i] == str2[i]:
            score.append(max(1 - abs(i-midpoint) ** 2 * 0.01, 0.05))
        else:
            score.append(0)
    summedscore = []
    for i in range(len1):
        badflag = 0
        subtotal = 0
        for j in range(i, len1):
            if score[j] != 0:
                subtotal += score[j]
            elif badflag > 0:
                summedscore.append(subtotal)
                break
            else:
                badflag += 1
        if badflag == 0:
            return subtotal
    return max(summedscore) if len(summedscore) > 0 else 0


def myscore4(str1, str2, brkpnt1, brkpnt2):
    if len(str1) * len(str2) == 0:
        return 0
    score = []
    brkpnt1 += 0.5
    brkpnt2 += 0.5
    for i in range(min(len(str1), len(str2))):
        if str1[i] == str2[i]:
            score.append(max(1 - max(abs(i-brkpnt1), abs(i-brkpnt2)) ** 1.5 * 0.01, 0.02))
        else:
            score.append(0)
    summedscore = []
    for i in range(len(str1)):
        badflag = 0
        subtotal = 0
        for j in range(i, len(str1)):
            if score[j] != 0:
                subtotal += score[j]
            elif badflag > 0:
                badflag += 1
                summedscore.append(subtotal)
                break
            else:
                badflag += 1
        if badflag <= 1:
            summedscore.append(subtotal)
    return max(summedscore) if len(summedscore) > 0 else 0
'''

def getline(thefilepath, desired_line_number):
    if desired_line_number < 1:
        return ''
    for current_line_number, line in enumerate(open(thefilepath, 'rU')):
        if current_line_number == desired_line_number - 1:
            return line
    return ''

'''
def EditDistance(word1, word2):
    m = len(word1) + 1
    n = len(word2) + 1
    dp = [[0 for i in range(n)] for j in range(m)]
    for i in range(n):
        dp[0][i] = i
    for i in range(m):
        dp[i][0] = i
    for i in range(1, m):
        for j in range(1, n):
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1,
                           dp[i - 1][j - 1] + (0 if word1[i - 1] == word2[j - 1] else 1))
    return dp[m - 1][n - 1]


def SWaligner(sequence1, sequence2):
    def get_match_score(s1, s2):
        if s1 == s2:
            return 3
        return -3

    def get_matrix_max(matrix):
        Max = matrix.max()
        for i in range(len(sequence2) + 1):
            for j in range(len(sequence1) + 1):
                if matrix[i][j] == Max:
                    return i, j

    s1 = ''
    s2 = ''
    gap = -2
    best_matrix = np.empty(shape=(len(sequence2) + 1, len(sequence1) + 1), dtype=int)
    for i in range(len(sequence2) + 1):
        for j in range(len(sequence1) + 1):
            if i == 0 or j == 0:
                best_matrix[i][j] = 0
            else:
                match = get_match_score(sequence2[i - 1], sequence1[j - 1])
                gap1_score = best_matrix[i - 1][j] + gap
                gap2_score = best_matrix[i][j - 1] + gap
                match_score = best_matrix[i - 1][j - 1] + match
                score = max(gap1_score, gap2_score, match_score)
                if score > 0:
                    best_matrix[i][j] = score
                else:
                    best_matrix[i][j] = 0
    # traceback
    i, j = get_matrix_max(best_matrix)
    besti, bestj = i, j
    while best_matrix[i][j] != 0:
        match = get_match_score(sequence2[i - 1], sequence1[j - 1])
        if i > 0 and j > 0 and best_matrix[i][j] == best_matrix[i - 1][j - 1] + match:
            s1 += sequence1[j - 1]
            s2 += sequence2[i - 1]
            i -= 1
            j -= 1
        elif i > 0 and best_matrix[i, j] == best_matrix[i - 1, j] + gap:
            s1 += '-'
            s2 += sequence2[i - 1]
            i -= 1
        else:
            s1 += sequence1[j - 1]
            s2 += '-'
            j -= 1
    return s1[::-1], s2[::-1], besti, bestj
'''

def findpos(chromo, pos):
    row = int(math.ceil(pos / reffilelinelength))
    col = pos % reffilelinelength
    if col == 0:
        col = reffilelinelength
    return row, col

'''
def getTCGAstring(chr, start, end):
    row1, col1 = findpos(chr, start)
    row2, col2 = findpos(chr, end)
    if row1 == row2:
        return ref[row1 - 1][col1 - 1:col2]
    resstring = ref[row1 - 1].rstrip()[col1 - 1:]
    for i in range(row1 + 1, row2):
        resstring += ref[i - 1].rstrip()
    return resstring + ref[row2 - 1][:col2]
'''


def getTCGAstring(chr, start, end):
    row1, col1 = findpos(chr, start)
    row2, col2 = findpos(chr, end)
    if row1 == row2:
        return refdict[chr][row1 - 1][col1 - 1:col2]
    resstring = refdict[chr][row1 - 1].rstrip()[col1 - 1:]
    for i in range(row1 + 1, row2):
        resstring += refdict[chr][i - 1].rstrip()
    return resstring + refdict[chr][row2 - 1][:col2]

'''
def getstring(row, col, halfwindowwidth, ref):
    try:
        if 1 + halfwindowwidth <= col <= 61 - halfwindowwidth:
            return ref[row-1][col-halfwindowwidth-1:col+halfwindowwidth-1]
        if col <= halfwindowwidth:
            return ref[row-2][59+col-halfwindowwidth:-1] + ref[row-1][:col+halfwindowwidth-1]
        return ref[row-1][col-halfwindowwidth-1:-1] + ref[row][:col+halfwindowwidth-61]
    except:
        print(str(row) + '\t' + str(col))
'''

def onethread(halfwindowwidth, fusionfile):
    geneset = {}
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
            genetype = ''
            if len(geneinfo) > 4:
                for item in geneinfo:
                    if item.find('gene_name') > -1:
                        genename = item[11:-1]
                    if item.find('gene_type') > -1:
                        genetype = item[11:-1]
                    if item.find('gene_biotype') > -1:
                        genetype = item[14:-1]
            else:
                continue
            if genetype == '' or genename == '':
                continue
            if genetype != 'protein_coding' and genetype != 'processed_transcript':
                continue
            if chromo not in geneset:
                geneset[chromo] = {}
            if genename not in geneset[chromo]:
                geneset[chromo][genename] = [[start, end], [[start, end]]]
            else:
                geneset[chromo][genename][1].append([start, end])
                geneset[chromo][genename][0][0] = min(geneset[chromo][genename][0][0], start)
                geneset[chromo][genename][0][1] = max(geneset[chromo][genename][0][1], end)
    refannotfile.close()
    allline = fusionfile.readlines()
    for line in allline:
        if line[0] == '#':
            continue
        try:
            info = line.split('\t')
            gene1 = info[0]
            gene2 = info[1]
            chr1 = info[3].replace(' ', '')
            chr2 = info[4].replace(' ', '')
            if not chr1.startswith('chr'):
                chr1 = 'chr' + chr1
            if not chr2.startswith('chr'):
                chr2 = 'chr' + chr2
            pos1 = int(info[5].replace(' ', ''))
            pos2 = int(info[6].replace(' ', ''))
            exoninterval1 = -1
            exoninterval2 = -1
            for g1 in geneset[chr1]:
                if geneset[chr1][g1][0][0] <= pos1 <= geneset[chr1][g1][0][1]:
                    gene1 = g1
                    break
            for g2 in geneset[chr2]:
                if geneset[chr2][g2][0][0] <= pos2 <= geneset[chr2][g2][0][1]:
                    gene2 = g2
                    break
            if gene1 not in geneset[chr1] or gene2 not in geneset[chr2]:
                pass
            else:
                for i in range(len(geneset[chr1][gene1][1])):
                    if geneset[chr1][gene1][1][i][0] <= pos1 <= geneset[chr1][gene1][1][i][1]:
                        exoninterval1 = i
                        break
                for i in range(len(geneset[chr2][gene2][1])):
                    if geneset[chr2][gene2][1][i][0] <= pos2 <= geneset[chr2][gene2][1][i][1]:
                        exoninterval2 = i
                        break
            # 获取外显子的序列
            if exoninterval2 == -1 or exoninterval1 == -1:
                res1d = getTCGAstring(chr1, pos1, pos1 + halfwindowslength)
                res1u = getTCGAstring(chr1, pos1 - halfwindowslength, pos1 - 1)
                res2d = getTCGAstring(chr2, pos2, pos2 + halfwindowslength)
                res2u = getTCGAstring(chr2, pos2 - halfwindowslength, pos2 - 1)
            else:
                res1d = ''
                res2d = ''
                res1u = ''
                res2u = ''
                last1 = 0
                last2 = 0
                if geneset[chr1][gene1][1][exoninterval1][1] - pos1 + 1 >= halfwindowslength + 1:
                    res1d = getTCGAstring(chr1, pos1, pos1 + halfwindowslength - 1)
                else:
                    res1d = getTCGAstring(chr1, pos1, geneset[chr1][gene1][1][exoninterval1][1])
                    last1 = geneset[chr1][gene1][1][exoninterval1][1]
                if geneset[chr2][gene2][1][exoninterval2][1] - pos1 + 1 >= halfwindowslength + 1:
                    res2d = getTCGAstring(chr2, pos2, pos2 + halfwindowslength - 1)
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
                    if geneset[chr1][gene1][1][fi1][1] - geneset[chr1][gene1][1][fi1][0] + 1 + len(
                            res1d) >= halfwindowslength + 1:
                        res1d += getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0],
                                               geneset[chr1][gene1][1][fi1][0] + halfwindowslength - len(res1d) - 1)
                    else:
                        res1d += getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0], geneset[chr1][gene1][1][fi1][1])
                        last1 = geneset[chr1][gene1][1][fi1][1]
                while len(res2d) < halfwindowslength:
                    fi2 += 1
                    if fi2 >= len(geneset[chr2][gene2][1]):
                        break
                    if last2 >= geneset[chr2][gene2][1][fi2][0]:
                        break
                    if geneset[chr2][gene2][1][fi2][1] - geneset[chr2][gene2][1][fi2][0] + 1 + len(
                            res2d) >= halfwindowslength + 1:
                        res2d += getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0],
                                               geneset[chr2][gene2][1][fi2][0] + halfwindowslength - len(res2d) - 1)
                    else:
                        res2d += getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0], geneset[chr2][gene2][1][fi2][1])
                        last2 = geneset[chr2][gene2][1][fi2][1]
                if - geneset[chr1][gene1][1][exoninterval1][0] + pos1 >= halfwindowslength:
                    res1u = getTCGAstring(chr1, pos1 - halfwindowslength, pos1 - 1)
                else:
                    res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][exoninterval1][0], pos1 - 1)
                    last1 = geneset[chr1][gene1][1][exoninterval1][0]
                if - geneset[chr2][gene2][1][exoninterval2][0] + pos2 >= halfwindowslength:
                    res2u = getTCGAstring(chr2, pos2 - halfwindowslength, pos2 - 1)
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
                    if geneset[chr1][gene1][1][fi1][1] - geneset[chr1][gene1][1][fi1][0] + 1 + len(
                            res1u) >= halfwindowslength:
                        res1u = getTCGAstring(chr1,
                                              geneset[chr1][gene1][1][fi1][1] - halfwindowslength + len(res1u) + 1,
                                              geneset[chr1][gene1][1][fi1][1]) + res1u
                    else:
                        res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0],
                                              geneset[chr1][gene1][1][fi1][1]) + res1u
                        last1 = geneset[chr1][gene1][1][fi1][0]
                while len(res2u) < halfwindowslength:
                    fi2 -= 1
                    if fi2 < 0:
                        break
                    if last2 <= geneset[chr2][gene2][1][fi2][1]:
                        break
                    if geneset[chr2][gene2][1][fi2][1] - geneset[chr2][gene2][1][fi2][0] + 1 + len(
                            res2u) >= halfwindowslength:
                        res2u = getTCGAstring(chr2,
                                              geneset[chr2][gene2][1][fi2][1] - halfwindowslength + len(res2u) + 1,
                                              geneset[chr2][gene2][1][fi2][1]) + res2u
                    else:
                        res2u = getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0],
                                              geneset[chr2][gene2][1][fi2][1]) + res2u
                        last2 = geneset[chr2][gene2][1][fi2][0]
            string1 = res1u.lower() + res1d.upper()
            string2 = res2u.lower() + res2d.upper()



            #string1 = getstring(row1, col1, halfwindowwidth, ref)
            #string2 = getstring(row2, col2, halfwindowwidth, ref)
            gccontent = GCcontent(string1, string2)
            '''
            string2c = ReverseComplement(string2)
            string1r = string1[::-1]
            SW1, SW2, bestj, besti = SWaligner(string1, string2)
            SW1c, SW2c, bestj2, besti2 = SWaligner(string1, string2c)
            SW1rc, SW2rc, bestj3, besti3 = SWaligner(string1r, string2c)
            SW1r, SW2r, bestj4, besti4 = SWaligner(string1r, string2)
            SW1t = [SW1, SW1c, SW1rc, SW1r]
            SW2t = [SW2, SW2c, SW2rc, SW2r]
            bestjt = [bestj, bestj2, bestj3, bestj4]
            bestit = [besti, besti2, besti3, besti4]
            stringlist1 = [string1, string1, string1r, string1r]
            stringlist2 = [string2, string2c, string2c, string2]
            brkpnt1 = []
            brkpnt2 = []
            for i in range(4):
                brkpnt1.append(Findbrkpnt_SWResult(SW1t[i], bestit[i], halfwindowwidth))
                brkpnt2.append(Findbrkpnt_SWResult(SW2t[i], bestjt[i], halfwindowwidth))
            # score = [myscore2(SW1, SW2), myscore2(SW1c, SW2c), myscore2(SW1rc, SW2rc), myscore2(SW1r, SW2r)]
            score = [myscore3(string1, string2), myscore3(string1, string2c), myscore3(string1r, string2c),
                    myscore3(string1r, string2)]
            # score = [myscore4(SW1, SW2, brkpnt1[0], brkpnt2[0]), myscore4(SW1c, SW2c, brkpnt1[1], brkpnt2[1]),
            #         myscore4(SW1rc, SW2rc, brkpnt1[2], brkpnt2[2]), myscore4(SW1r, SW2r, brkpnt1[3], brkpnt2[3])]
            bigindex = score.index(max(score))
            if SW1t[bigindex] == '':
                SW1t[bigindex] = '---'
            if SW2t[bigindex] == '':
                SW2t[bigindex] = '---'
            print(line[:-1] + '\t' + SW1t[bigindex] + '\t' + SW2t[bigindex] + '\t' + str(max(score)) + '\t' + str(gccontent[0])
                  + '\t' + str(gccontent[1]) + '\t' + str(brkpnt1[bigindex]) + '\t' + str(brkpnt2[bigindex]) + '\t' +
                  stringlist1[bigindex] + '\t' + stringlist2[bigindex])
            '''
            print(line.rstrip() + '\t-\t-\t0\t' + str(gccontent[0]) + '\t' + str(gccontent[1]) + '\t0\t0\t' + string1 + '\t' + string2)
        except:
            sys.stderr.write('Bad Line:' + line)
'''
def Simulation_random(length):
    string1 = simulateString(length)
    string2 = simulateString(length)
    string2c = ReverseComplement(string2)
    SW1, SW2, bestj, besti = SWaligner(string1, string2)
    SW1c, SW2c, bestj2, besti2 = SWaligner(string1, string2c)
    score = myscore2(SW1, SW2)
    score2 = myscore2(SW1c, SW2c)
    if SW1 == '':
        SW1 = '---'
    if SW2 == '':
        SW2 = '---'
    if SW1c == '':
        SW1c = '---'
    if SW2c == '':
        SW2c = '---'
    if score >= score2:
        brkpnt1 = Findbrkpnt_SWResult(SW1, besti, int(length / 2))
        brkpnt2 = Findbrkpnt_SWResult(SW2, bestj, int(length / 2))
        print(string1 + '\t' + string2 + '\t' + SW1 + '\t' + SW2 + '\t' + str(score) + '\t' + str(brkpnt1) + '\t' + str(brkpnt2))
    else:
        brkpnt1 = Findbrkpnt_SWResult(SW1c, besti2, int(length / 2))
        brkpnt2 = Findbrkpnt_SWResult(SW2c, bestj2, int(length / 2))
        print(string1 + '\t' + string2c + '\t' + SW1c + '\t' + SW2c + '\t' + str(score2) + '\t' + str(brkpnt1) + '\t' + str(brkpnt2))


def Simulation_trans(length, ref):
    lines = random.sample(ref, 10000)
    line1 = lines[:5000]
    line2 = lines[5000:]
    for i in range(5000):
        if line1[i][0] == '>' or line2[i][0] == '>':
            continue
        if line1[i].find('N') > -0.5 or line2[i].find('N') > -0.5:
            continue
        len1 = len(line1[i][:-1])
        len2 = len(line2[i][:-1])
        if len1 <= length or len2 <= length:
            continue
        pos1 = random.randint(0, len1-length)
        string1 = line1[i][pos1:pos1+length]
        pos2 = random.randint(0, len2 - length)
        string2 = line2[i][pos2:pos2+length]
        string2c = ReverseComplement(string2)
        SW1, SW2, bestj, besti = SWaligner(string1, string2)
        SW1c, SW2c, bestj2, besti2 = SWaligner(string1, string2c)
        score = myscore2(SW1, SW2)
        score2 = myscore2(SW1c, SW2c)
        if SW1 == '':
            SW1 = '---'
        if SW2 == '':
            SW2 = '---'
        if SW1c == '':
            SW1c = '---'
        if SW2c == '':
            SW2c = '---'
        if score >= score2:
            brkpnt1 = Findbrkpnt_SWResult(SW1, besti, int(length / 2))
            brkpnt2 = Findbrkpnt_SWResult(SW2, bestj, int(length / 2))
            print(string1 + '\t' + string2 + '\t' + SW1 + '\t' + SW2 + '\t' + str(score) + '\t' + str(
                brkpnt1) + '\t' + str(brkpnt2))
        else:
            brkpnt1 = Findbrkpnt_SWResult(SW1c, besti2, int(length / 2))
            brkpnt2 = Findbrkpnt_SWResult(SW2c, bestj2, int(length / 2))
            print(string1 + '\t' + string2c + '\t' + SW1c + '\t' + SW2c + '\t' + str(score2) + '\t' + str(
                brkpnt1) + '\t' + str(brkpnt2))
'''

def main():
    onethread(halfwindowwidth, fusionfile)


if __name__ == "__main__":
    #rownum = [1, 4154180, 8207504, 11507879, 14693785, 17709041, 20560960, 23213273, 25652675, 28006234, 30265148, 32515258,34746124, 36665623, 38454783, 40163641, 41669555, 43022810, 44324099, 45309584, 46360011, 47162177, 48017255, 50605099]
    halfwindowwidth = 10
    halfwindowslength = 100
    fusionfile = open(sys.argv[1])
    reffilepath = sys.argv[2]
    reference = open(reffilepath)
    refannotfile = open(sys.argv[3])
    #ref = reference.readlines()
    currentchr = ''
    refdict = {}
    for line in reference.readlines():
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
    main()
