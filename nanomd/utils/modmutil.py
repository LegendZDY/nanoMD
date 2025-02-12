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
import argparse
import re
import gzip
import pysam
import concurrent.futures

class getModifications:
    BASE_PATTERNS = {
        'A': re.compile(r'A'),
        'T': re.compile(r'T'),
        'C': re.compile(r'C')
    }

    def __init__(self, inputFq, samFile, bedFile, output, workers=4):
        self.inputFq = inputFq
        self.samFile = samFile
        self.bedFile = bedFile
        self.output = output
        self.workers = workers

    def getAnnotation(self):
        def read_bed_file():
            with open(self.bedFile, 'r') as f:
                for line in f:
                    line = line.strip().split("\t")
                    gid = line[0] + ":" + line[1]
                    enst = line[3]
                    base = line[4]
                    yield gid, (enst, base)
        return dict(read_bed_file())

    def posFind(self, modList):
        counter = 0
        posList = []
        for mod in modList:
            counter += int(mod)
            posList.append(counter)
            counter += 1
        return posList

    # 找到修饰位置和对应的p值
    def intersectionPosition(self, mm, baseList, mlNum):
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
        return readPos, readPosPvalue


    def process_record(self, name, seq):
        abase = re.compile(r'A')
        tbase = re.compile(r'T')
        cbase = re.compile(r'C')
        namelist = name.split("\t")
        id = namelist[0].split("@")[-1]  # 获取ID
        MMtagList = namelist[-2].split("MM:Z:")[-1].split(";")
        MLtag = namelist[-1].split(",")[1:]
        Alist = [match.start() + 1 for match in abase.finditer(seq)]
        Tlist = [match.start() + 1 for match in tbase.finditer(seq)]
        Clist = [match.start() + 1 for match in cbase.finditer(seq)]
        mlstart = 0
        rpos_mod = {}
        for mm in MMtagList:
            mmlen = len(mm.split(","))
            if mmlen > 1:
                mlend = mlstart + mmlen - 1
                mlNum = MLtag[mlstart:mlend]
                mlstart += mmlen - 1
                mod_map = {
                    'A+a': Alist,
                    'A+17596': Alist,
                    'T+17802': Tlist,
                    'C+m': Clist
                }
                for mod, base_list in mod_map.items():
                    if mm.startswith(mod):
                        readPos, readPosPvalue = self.intersectionPosition(mm, base_list, mlNum)
                        for rpos, pvalue in zip(readPos, readPosPvalue):
                            if pvalue > 0.98:
                                if rpos_mod.get(id) is None:
                                    rpos_mod[id] = [(rpos, pvalue, mod)]
                                else:
                                    rpos_mod[id].append((rpos, pvalue, mod))
        return rpos_mod

    def get_mod_position_with_read(self):
        rpos_mod = {}
        with gzip.open(self.inputFq, 'rt') as f:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                futures = []
                while True:
                    try:
                        # 读取四行fastq格式数据
                        name = next(f).strip()
                        seq = next(f).strip()
                        plus = next(f).strip()
                        qual = next(f).strip()
                        futures.append(executor.submit(self.process_record, name, seq))
                    except StopIteration:
                        break

                for future in concurrent.futures.as_completed(futures):
                    result = future.result()
                    for id, mod_info in result.items():
                        if id in rpos_mod:
                            rpos_mod[id].extend(mod_info)
                        else:
                            rpos_mod[id] = mod_info

        return rpos_mod


    def get_mod_position_with_sam(self):
        rpos_mod = self.get_mod_position_with_read()
        annot = self.getAnnotation()
        with open(self.output, "w") as f:
            for read in pysam.AlignmentFile(self.samFile, "r"):
                if not read.is_unmapped and rpos_mod.get(read.query_name) is not None:
                    pos = read.reference_start
                    read_pos = 1
                    strand = "-" if read.is_reverse else "+"
                    # 循环遍历 cigartuples，将读取的读和基因组位置合并后直接写入文件
                    for op, length in read.cigartuples:
                        if op == 0:  # M: 匹配或错配
                            # 匹配的碱基在基因组上的位置
                            for i in range(length):
                                mod_list = rpos_mod[read.query_name]
                                for rpos, pvalue, mod in mod_list:
                                    if rpos == read_pos + i:
                                        genome_id = f"chr{read.reference_name}:{pos + i}"
                                        if annot.get(genome_id) is not None:
                                            enst, base = annot[genome_id]
                                            f.write(f"{read.reference_name}\t{pos + i}\t{pos + i+1}\t{pvalue}\t{enst}\t{strand}\t{mod}\t{base}\n")
                            pos += length
                            read_pos += length
                        elif op == 1:  # I: 插入
                            # 插入的碱基在参考基因组上不占位置，但在 read 上有位置
                            read_pos += length
                        elif op == 2:  # D: 删除
                            # 删除的碱基在参考基因组上仍占位置，但跳过它们
                            pos += length
                        elif op == 4:  # S: 跳过
                            # 跳过的碱基不占基因组位置，但在 read 上有位置
                            read_pos += length
                        elif op == 5:  # H: 硬剪切
                            # 硬剪切不占基因组位置
                            read_pos += length
                        else:
                            print(f"op: {op}, length: {length} not supported")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parallel modifacts detection from raw fastq data.'
    )
    parser.add_argument('--version', action='version', version='modkit 1.0')
    parser.add_argument('-i', '--inputFq', required=True)
    parser.add_argument('-s', '--samFile', required=True)
    parser.add_argument('-b', '--bedFile', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-w', '--workers', default=4, type=int)
    args = parser.parse_args()
    mod = getModifications(args.inputFq, args.samFile, args.bedFile, args.output, args.workers)
    mod.get_mod_position_with_sam()
