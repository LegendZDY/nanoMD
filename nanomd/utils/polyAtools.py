# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250612
Description: This file contains the functions for converting pod5 files to fast5 files and indexing the fastq file using nanopolish index.
"""
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List
from more_itertools import chunked

import pysam
import pod5 as p5
from pod5.tools.utils import DEFAULT_THREADS, collect_inputs, limit_threads
from pod5.tools.pod5_convert_to_fast5 import *
from basebio import run_command

def convert_to_fast5_with_summary_file(
    inputs: List[Path],
    output: Path,
    summary_file: str,
    fastq_name: str,
    recursive: bool = False,
    threads: int = DEFAULT_THREADS,
    force_overwrite: bool = False,
    file_read_count: int = 4000,
):
    if output.exists() and not output.is_dir():
        raise FileExistsError("Cannot output to a file")

    threads = limit_threads(threads)
    columns = ["filename_fastq", "filename_fast5", "read_id", "run_id", "channel", "mux", "start_time", "duration", "num_events", "passes_filtering", "template_start", "num_events_template", "template_duration", "sequence_length_template", "mean_qscore_template", "strand_score_template", "median_template", "mad_template", "pore_type", "experiment_id", "sample_id", "end_reason", "alias", "type", "barcode_arrangement", "barcode_full_arrangement", "barcode_kit", "barcode_variant", "barcode_score", "barcode_front_id", "barcode_front_score", "barcode_front_refseq", "barcode_front_foundseq", "barcode_front_foundseq_length", "barcode_front_begin_index", "barcode_rear_id", "barcode_rear_score", "barcode_rear_refseq", "barcode_rear_foundseq", "barcode_rear_foundseq_length", "barcode_rear_end_index"]

    # 数据行（单行列表）
    data_row = ["test", "885", "1", "56.132750", "0.756500", "605", "TRUE", "56.154000", "588", "0.735250", "153", "8.314826", "0.000000", "94.303551", "8.406906", "not_set", "20220426-CDS0204-PAK00515", "no_sample", "signal_positive", "barcode11", "test_sample", "barcode11", "LWB11_var1", "LWB", "var1", "92.750000", "BC11_FWD", "89.250000", "CCGTGACGTT", "51", "34", "BC11_REV", "92.750000", "GCAATAA", "GCACATC", "52", "20"]
    
    with open(summary_file, "w") as fo:
        fo.write("\t".join(columns) + "\n")
        with ProcessPoolExecutor(max_workers=threads) as executor:
            total_reads = 0
            futures: Dict[Future, Path] = {}
            # Enumerate over input pod5 files
            for input_idx, source in enumerate(
                collect_inputs(inputs, recursive, "*.pod5", threads=threads)
            ):
                # Open the inputs to read the read ids
                with p5.Reader(source) as reader:
                    for chunk_idx, read_ids in enumerate(
                        chunked(reader.read_ids, file_read_count)
                    ):
                        dest = (
                            output / f"{source.stem}.{chunk_idx}_{input_idx}.fast5"
                        ).resolve()

                        if dest.exists() and not force_overwrite:
                            raise FileExistsError(
                                "Output path points to an existing file and --force-overwrite not set"
                            )

                        kwargs = {
                            "source": source,
                            "dest": dest,
                            "read_ids": read_ids,
                        }
                        futures[executor.submit(convert_pod5_to_fast5, **kwargs)] = dest  # type: ignore

                        for read_id in read_ids:
                            raw = [fastq_name, f"{source.stem}.{chunk_idx}_{input_idx}.fast5", read_id] + data_row
                            fo.write("\t".join(raw) + "\n")

                    total_reads += len(reader.read_ids)

        print(f"Converting pod5s into {len(futures)} fast5 files. Please wait...")

        status = StatusMonitor(file_count=len(inputs))
        status.increment(files_started=len(inputs), read_count=total_reads)

        for idx, future in enumerate(as_completed(futures)):
            (reads_converted, samples_converted) = future.result()

            status.increment(
                files_ended=1,
                sample_count=samples_converted,
                reads_processed=reads_converted,
            )
            status.print_status()

        status.print_status(force=True)

        print("Conversion complete")

def index_fastq(fast5_dir, summary_file, fastq):
    """
    Index the fastq file using nanopolish index.
    
    Args:
        fast5_dir (str): The directory where the fast5 files are located.
        summary_file (str): The path to the sequencing summary file.
        fastq (str): The path to the fastq file.
    
    Examples:
        index_fastq("fast5_dir", "summary.txt", "fastq.gz")
    """

    cmd = f"nanopolish index --directory={fast5_dir} --sequencing-summary={summary_file} {fastq}".split()
    run_command(cmd)

def detect_polyA(fastq, bam, transcriptome, output_ploya, threads=8):
    """
    Detect polyA using nanopolish polya.

    Args:
        fastq (str): The path to the fastq file.
        bam (str): The path to the bam file.
        transcriptome (str): The path to the transcriptome file.
        output_ploya (str): The path to the output polyA file.
        threads (int): The number of threads to use.

    Examples:
        detect_polyA("fastq.gz", "bam", "transcriptome.fa", "polyA.tsv", threads=8)
    """

    cmd = f"nanopolish polya --threads={threads} --reads={fastq} --bam={bam} --genome={transcriptome} > {output_ploya}"
    run_command(cmd, use_shell=True)

class PolyADetector:
    """
    polyA detector.

    Args:
        bam_path (str): The path to the bam file.
        output_path (str): The path to the output file.
        min_a_length (int): The minimum length of polyA.
        max_non_a (int): The maximum number of non-A bases allowed in polyA.
    
    Examples:
        detector = PolyADetector("bam", "polyA.tsv", min_a_length=10, max_non_a=3)
        detector.analyze()
    """
    def __init__(self, bam_path, output_path, min_a_length=6, max_non_a=3):
        self.bam_path = bam_path
        self.output_path = output_path
        self.min_a_length = min_a_length
        self.max_non_a = max_non_a
    
    def find_longest_polyA(self, seq):
        """
        Use sliding window algorithm to find the longest consecutive A sequence, allowing a few non-A bases.
        """
        if not seq:
            return "", 0, 0, 0.0
        
        n = len(seq)
        left = 0
        non_a_count = 0
        max_length = 0
        a_count_in_max = 0
        best_start = 0
        best_end = -1
        
        for right in range(n):
            # 更新非A碱基计数
            if seq[right].upper() != 'A':
                non_a_count += 1
            
            # 移动左指针直到非A碱基数量在允许范围内
            while non_a_count > self.max_non_a:
                if seq[left].upper() != 'A':
                    non_a_count -= 1
                left += 1
            
            # 计算当前窗口长度和A碱基数量
            current_length = right - left + 1
            current_a_count = current_length - non_a_count
            
            # 检查当前窗口是否更长
            if current_length > max_length:
                max_length = current_length
                a_count_in_max = current_a_count
                best_start = left
                best_end = right
        
        # 计算A的比例
        a_ratio = a_count_in_max / max_length if max_length > 0 else 0.0
        
        # 提取polyA序列
        polyA_seq = seq[best_start:best_end+1] if max_length > 0 else ""
        
        return polyA_seq, max_length, a_count_in_max, a_ratio
    
    def process_read(self, read):
        # 跳过未比对的read
        if read.is_unmapped or not read.cigartuples:
            return None
        
        read_name = read.query_name
        seq = read.query_sequence
        if not seq:  # 跳过空序列
            return None
            
        flag = read.flag
        cigar = read.cigartuples
        
        # 确定需要检查的CIGAR操作位置
        if flag == 0:    # 正向比对
            target_op = cigar[-1]  # 最后一个操作
        elif flag == 16: # 反向比对
            target_op = cigar[0]   # 第一个操作
        else:
            return None  # 跳过其他比对方向
        
        op_type, op_len = target_op
        
        # 只处理软裁剪操作（S=4）
        if op_type != 4:
            return None
        
        # 提取待检测序列
        if flag == 0:
            # 正向比对：取序列末尾的软裁剪区域
            target_seq = seq[-op_len:] if op_len <= len(seq) else seq
        else:  # flag == 16
            # 反向比对：取序列开头的软裁剪区域
            target_seq = seq[:op_len] if op_len <= len(seq) else seq
        
        # 在目标序列中查找最长连续A（允许少量非A）
        polyA_seq, total_length, a_count, a_ratio = self.find_longest_polyA(target_seq)
        
        # 判断是否存在polyA
        has_polyA = a_count >= self.min_a_length
        
        return {
            "read_name": read_name,
            "ref_name": read.reference_name,
            "strand": "+" if flag == 0 else "-",
            "polyA_seq": polyA_seq,
            "polyA_region_length": total_length,
            "a_count": a_count,
            "a_ratio": a_ratio,
            "has_polyA": "Yes" if has_polyA else "No"
        }
    
    def analyze(self):
        with open(self.output_path, 'w') as fout:
            fout.write("readName\trefName\tstrand\tpolyASeq\tpolyALength\tACount\tARatio\tHasPolyA\n")
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                for read in bam:
                    result = self.process_read(read)
                    if result:
                        fout.write("\t".join([
                            result["read_name"],
                            result["ref_name"],
                            result["strand"],
                            result["polyA_seq"],
                            str(result["polyA_region_length"]),
                            str(result["a_count"]),
                            f"{result['a_ratio']:.3f}",
                            result["has_polyA"]
                        ]) + "\n")