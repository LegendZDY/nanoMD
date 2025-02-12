#!/usr/bin/python3
# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250210
Description:
This Script of raw fq data from Legendzdy.  Of course, this means there's a possibility
for other ways.Use at your own discretion.
"""
#%%fowardDir
import argparse, re, gzip
import pysam

class getModifications(object):
    def __init__(self, inputFq, samFile, bedFile, output):
        self.inputFq = inputFq
        self.samFile = samFile
        self.bedFile = bedFile
        self.output = output
    
    def posFind(self, modList):
        counter = 0
        posList = []
        for mod in modList:
            counter += int(mod)
            posList.append(counter)
            counter += 1
        return posList
    
    def getAnnotation(self):
        annot = {}
        with open(self.bedFile, 'r') as f:
            for line in f:
                line = line.strip().split("\t")
                gid = line[0] + ":" + line[1]
                enst = line[3]
                base = line[4]
                annot[gid] = [enst, base]
        return annot
    
    def getModPositionWithReadAndGenome(self):
        samfile = pysam.AlignmentFile(self.samFile, "r")
        readAndGenome = {}
        for read in samfile:
            if not read.is_unmapped:
                genome_positions = []
                read_positions = []
                pos = read.reference_start
                read_pos = 1
                for op, length in read.cigartuples:
                    # print(f"op: {op}, length: {length}")
                    if op == 0:  # M: 匹配或错配
                        # 匹配的碱基在基因组上的位置
                        genome_positions.extend(range(pos, pos + length))
                        read_positions.extend(range(read_pos, read_pos + length))
                        pos += length
                        read_pos += length
                    elif op == 1:  # I: 插入
                        # 插入的碱基在参考基因组上不占位置，但在 read 上有位置
                        # read_positions.extend(range(read_pos, read_pos + length))
                        read_pos += length
                    elif op == 2:  # D: 删除
                        # 删除的碱基在参考基因组上仍占位置，但跳过它们
                        # genome_positions.extend(range(pos, pos + length))
                        pos += length
                    elif op == 4:  # S: 跳过
                        # 跳过的碱基不占基因组位置，但在 read 上有位置
                        # read_positions.extend(range(read_pos, read_pos + length))
                        read_pos += length
                    elif op == 5:  # H: 硬剪切
                        # 硬剪切不占基因组位置
                        read_pos += length
                    else:
                        print(f"op: {op}, length: {length} not supported")
                if read.is_reverse:
                    strand = "-"
                else:
                    strand = "+"
                readAndGenome[read.query_name] = [read_positions, genome_positions, read.reference_name, strand]
        return readAndGenome
    
    def intersectionPosition(self, mm, baseList, mlNum, refReadPos, refGenomePos):
        modList = mm.split(",")[1:]
        posList = self.posFind(modList)
        mlcount = 0
        readPos = []
        readPosPvalue = []
        for pos in posList:
            rpos = baseList[pos]
            pvalue = round(int(mlNum[mlcount])/255, 3)
            readPos.append(rpos)
            readPosPvalue.append(pvalue)
            mlcount += 1
        intersection = list(set(refReadPos) & set(readPos))
        if len(intersection) != 0:
            PosRead = [readPos.index(x) for x in intersection]
            PosGeno = [refReadPos.index(x) for x in intersection]
            targetPvalue = [readPosPvalue[i] for i in PosRead]
            targetGeno = [refGenomePos[i] for i in PosGeno]
        else:
            targetGeno = []
            targetPvalue = []
        return targetGeno, targetPvalue
    
    def getModPositionWithRead(self):
        annot = self.getAnnotation()
        readAndGenome = self.getModPositionWithReadAndGenome()
        with open(self.output, 'w') as out:
            abase = re.compile(r'A')
            tbase = re.compile(r'T')
            cbase = re.compile(r'C')
            with gzip.open(self.inputFq, 'rt') as f:
                while True:
                    try:
                        name = next(f).strip()
                        seq = next(f).strip()
                        plus = next(f).strip()
                        qual = next(f).strip()
                        
                        namelist = name.split("\t")
                        id = namelist[0].split("@")[-1]
                        if id in readAndGenome:
                            refReadPos = readAndGenome[id][0]
                            refGenomePos = readAndGenome[id][1]
                            chrom = readAndGenome[id][2]
                            strand = readAndGenome[id][3]
                            MMtagList = namelist[-2].split("MM:Z:")[-1].split(";")
                            # print(MMtagList)
                            MLtag = namelist[-1].split(",")[1:]
                            Alist = [match.start() + 1 for match in abase.finditer(seq)]
                            Tlist = [match.start() + 1 for match in tbase.finditer(seq)]
                            Clist = [match.start() + 1 for match in cbase.finditer(seq)]
                            mlstart = 0
                            for mm in MMtagList:
                                mmlen = len(mm.split(","))
                                if mmlen > 1:
                                    mlend = mlstart + mmlen - 1
                                    mlNum = MLtag[mlstart:mlend]
                                    mlstart += mmlen -1
                                    if mm.startswith('A+a'):
                                        targetGeno, targetPvalue = self.intersectionPosition(mm, Alist, mlNum, refReadPos, refGenomePos)
                                        if len(targetGeno) != 0:
                                            for i in range(len(targetGeno)):
                                                start = targetGeno[i]
                                                end = start + 1
                                                gid = "chr" + chrom + ":" + str(start)
                                                modType = "m6A"
                                                if gid in annot:
                                                    enst = annot[gid][0]
                                                    base = annot[gid][1]
                                                    out.write(f"chr{chrom}\t{start}\t{end}\t{targetPvalue[i]}\t{enst}\t{strand}\t{modType}\t{base}\n")
                                                else:
                                                    print(f"gid: {gid} not in annot")
                                    elif mm.startswith('A+17596') :
                                        targetGeno, targetPvalue = self.intersectionPosition(mm, Alist, mlNum, refReadPos, refGenomePos)
                                        if len(targetGeno) != 0:
                                            for i in range(len(targetGeno)):
                                                start = targetGeno[i]
                                                end = start + 1
                                                gid = "chr" + chrom + ":" + str(start)
                                                modType = "m1A"
                                                if gid in annot:
                                                    enst = annot[gid][0]
                                                    base = annot[gid][1]
                                                    out.write(f"chr{chrom}\t{start}\t{end}\t{targetPvalue[i]}\t{enst}\t{strand}\t{modType}\t{base}\n")
                                                else:
                                                    print(f"gid: {gid} not in annot")
                                    elif mm.startswith('T+17802'):
                                        targetGeno, targetPvalue = self.intersectionPosition(mm, Tlist, mlNum, refReadPos, refGenomePos)
                                        if len(targetGeno) != 0:
                                            for i in range(len(targetGeno)):
                                                start = targetGeno[i]
                                                end = start + 1
                                                gid = "chr" + chrom + ":" + str(start)
                                                modType = "psi"
                                                if gid in annot:
                                                    enst = annot[gid][0]
                                                    base = annot[gid][1]
                                                    out.write(f"chr{chrom}\t{start}\t{end}\t{targetPvalue[i]}\t{enst}\t{strand}\t{modType}\t{base}\n")
                                                else:
                                                    print(f"gid: {gid} not in annot")
                                    elif mm.startswith('C+m'):
                                        targetGeno, targetPvalue = self.intersectionPosition(mm, Clist, mlNum, refReadPos, refGenomePos)
                                        if len(targetGeno) != 0:
                                            for i in range(len(targetGeno)):
                                                start = targetGeno[i]
                                                end = start + 1
                                                gid = "chr" + chrom + ":" + str(start)
                                                modType = "m5C"
                                                if gid in annot:
                                                    enst = annot[gid][0]
                                                    base = annot[gid][1]
                                                    out.write(f"chr{chrom}\t{start}\t{end}\t{targetPvalue[i]}\t{enst}\t{strand}\t{modType}\t{base}\n")
                                                else:
                                                    print(f"gid: {gid} not in annot")         
                    except StopIteration:
                        break
                    
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='modifacts detection from raw fastq data.'
    )
    parser.add_argument('--version', action='version', version='modkit Version: %(prog)s 0.23')
    parser.add_argument('-i', '--inputFq',  required=True, help='The input fastq file with MM tag.')
    parser.add_argument('-s', '--samFile',  required=True, help='The SAM file formed by mapping the input fastq file.')
    parser.add_argument('-b', '--bedFile',  required=True, help='The annotated BED file.')
    parser.add_argument('-o', '--output',  required=True, help='The output file.')
    args = parser.parse_args()
    mod = getModifications(args.inputFq, args.samFile, args.bedFile, args.output)
    mod.getModPositionWithRead()