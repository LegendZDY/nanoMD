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

import argparse, re, gzip
import pysam
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import sqlite3
import os

class getModifications(object):
    def __init__(self, inputFq, samFile, bedFile, output, max_processes):
        self.inputFq = inputFq
        self.samFile = samFile
        self.bedFile = bedFile
        self.output = output
        self.max_processes = max_processes
        self.max_threads = max_processes
        self.db_file = 'annotations.db'
        self.create_database()
    
    def create_database(self):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS annotations
                     (gid TEXT PRIMARY KEY, enst TEXT, base TEXT)''')
        c.execute('''CREATE TABLE IF NOT EXISTS readAndGenome
                     (query_name TEXT PRIMARY KEY, read_positions TEXT, genome_positions TEXT, reference_name TEXT, strand TEXT)''')
        conn.commit()
        conn.close()
        self.load_annotations()
    
    def load_annotations(self):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        with open(self.bedFile, 'r') as f:
            for line in f:
                line = line.strip().split("\t")
                gid = line[0] + ":" + line[1]
                enst = line[3]
                base = line[4]
                c.execute("INSERT OR REPLACE INTO annotations (gid, enst, base) VALUES (?, ?, ?)", (gid, enst, base))
        conn.commit()
        conn.close()
    
    def get_annotation(self, gid):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        c.execute("SELECT enst, base FROM annotations WHERE gid=?", (gid,))
        result = c.fetchone()
        conn.close()
        return result
    
    def get_read_and_genome(self, query_name):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        c.execute("SELECT read_positions, genome_positions, reference_name, strand FROM readAndGenome WHERE query_name=?", (query_name,))
        result = c.fetchone()
        conn.close()
        return result
    
    def posFind(self, modList):
        counter = 0
        posList = []
        for mod in modList:
            counter += int(mod)
            posList.append(counter)
            counter += 1
        return posList
    
    def getModPositionWithReadAndGenome(self):
        # 使用ProcessPoolExecutor进行并行处理
        with ProcessPoolExecutor(max_workers=self.max_processes) as executor:
            futures = []
            for read in pysam.AlignmentFile(self.samFile, "r"):
                if not read.is_unmapped:
                    future = executor.submit(self.process_read, read)
                    futures.append(future)
            
            results = [future.result() for future in futures]
        
        # 将结果写入数据库
        self.save_read_and_genome(results)
    
    def process_read(self, read):
        genome_positions = []
        read_positions = []
        pos = read.reference_start
        read_pos = 1
        for op, length in read.cigartuples:
            if op == 0:  # M: 匹配或错配
                genome_positions.extend(range(pos, pos + length))
                read_positions.extend(range(read_pos, read_pos + length))
                pos += length
                read_pos += length
            elif op == 1:  # I: 插入
                read_pos += length
            elif op == 2:  # D: 删除
                pos += length
            elif op == 4:  # S: 跳过
                read_pos += length
            elif op == 5:  # H: 硬剪切
                read_pos += length
            else:
                print(f"unknown cigar operation: {op}")
        
        strand = "-" if read.is_reverse else "+"
        
        return {
            'query_name': read.query_name,
            'read_positions': ','.join(map(str, read_positions)),
            'genome_positions': ','.join(map(str, genome_positions)),
            'reference_name': read.reference_name,
            'strand': strand
        }
    
    def save_read_and_genome(self, results):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        c.executemany("INSERT OR REPLACE INTO readAndGenome (query_name, read_positions, genome_positions, reference_name, strand) VALUES (?, ?, ?, ?, ?)",
                      [(result['query_name'], result['read_positions'], result['genome_positions'], result['reference_name'], result['strand']) for result in results])
        conn.commit()
        conn.close()
    
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
    
    def process_fastq_entry(self, entry):
        name, seq, plus, qual = entry
        namelist = name.split("\t")
        id = namelist[0].split("@")[-1]
        readAndGenome = self.get_read_and_genome(id)
        if readAndGenome:
            refReadPos = list(map(int, readAndGenome[0].split(',')))
            refGenomePos = list(map(int, readAndGenome[1].split(',')))
            chrom = readAndGenome[2]
            strand = readAndGenome[3]
            MMtagList = namelist[-2].split("MM:Z:")[-1].split(";")
            MLtag = namelist[-1].split(",")[1:]
            Alist = [match.start() + 1 for match in re.finditer('A', seq)]
            Tlist = [match.start() + 1 for match in re.finditer('T', seq)]
            Clist = [match.start() + 1 for match in re.finditer('C', seq)]
            mlstart = 0
            results = []
            for mm in MMtagList:
                mmlen = len(mm.split(","))
                if mmlen > 1:
                    mlend = mlstart + mmlen - 1
                    mlNum = MLtag[mlstart:mlend]
                    mlstart += mmlen - 1
                    if mm.startswith('A+a'):
                        targetGeno, targetPvalue = self.intersectionPosition(mm, Alist, mlNum, refReadPos, refGenomePos)
                        results.extend(self.generate_output(targetGeno, targetPvalue, chrom, strand, 'm6A'))
                    elif mm.startswith('A+17596'):
                        targetGeno, targetPvalue = self.intersectionPosition(mm, Alist, mlNum, refReadPos, refGenomePos)
                        results.extend(self.generate_output(targetGeno, targetPvalue, chrom, strand, 'm1A'))
                    elif mm.startswith('T+17802'):
                        targetGeno, targetPvalue = self.intersectionPosition(mm, Tlist, mlNum, refReadPos, refGenomePos)
                        results.extend(self.generate_output(targetGeno, targetPvalue, chrom, strand, 'psi'))
                    elif mm.startswith('C+m'):
                        targetGeno, targetPvalue = self.intersectionPosition(mm, Clist, mlNum, refReadPos, refGenomePos)
                        results.extend(self.generate_output(targetGeno, targetPvalue, chrom, strand, 'm5C'))
            return results
        return []

    def generate_output(self, targetGeno, targetPvalue, chrom, strand, modType):
        outputs = []
        for i in range(len(targetGeno)):
            start = targetGeno[i]
            end = start + 1
            gid = f"chr{chrom}:{start}"
            annot = self.get_annotation(gid)
            if annot:
                enst, base = annot
                outputs.append(f"chr{chrom}\t{start}\t{end}\t{targetPvalue[i]}\t{enst}\t{strand}\t{modType}\t{base}\n")
            else:
                print(f"gid: {gid} not in annot")
        return outputs
    
    def getModPositionWithRead(self):
        self.getModPositionWithReadAndGenome()
        
        with open(self.output, 'w') as out, ThreadPoolExecutor(max_workers=self.max_threads) as executor:
            futures = []
            with gzip.open(self.inputFq, 'rt') as f:
                while True:
                    try:
                        name = next(f).strip()
                        seq = next(f).strip()
                        plus = next(f).strip()
                        qual = next(f).strip()
                        entry = (name, seq, plus, qual)
                        future = executor.submit(self.process_fastq_entry, entry)
                        futures.append(future)
                    except StopIteration:
                        break
            
            for future in futures:
                results = future.result()
                for line in results:
                    out.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='modifacts detection from raw fastq data.'
    )
    parser.add_argument('--version', action='version', version='modkit Version: %(prog)s 0.23')
    parser.add_argument('-i', '--inputFq', required=True, help='The input fastq file with MM tag.')
    parser.add_argument('-s', '--samFile', required=True, help='The SAM file formed by mapping the input fastq file.')
    parser.add_argument('-b', '--bedFile', required=True, help='The annotated BED file.')
    parser.add_argument('-o', '--output', required=True, help='The output file.')
    parser.add_argument('-p', '--max-processes', type=int, default=4, help='Maximum number of processes for processing SAM file (default: 4).')
    args = parser.parse_args()
    mod = getModifications(args.inputFq, args.samFile, args.bedFile, args.output, args.max_processes)
    mod.getModPositionWithRead()