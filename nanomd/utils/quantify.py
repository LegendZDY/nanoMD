# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250612
Description: This file contains the quantification function.
"""
import os, glob
import pandas as pd
from collections import OrderedDict
from basebio import run_command, minimap2

def salmon_map(input, reference, output):
    """
    Map with minimap2.
    Args:
        input: input fasta file.
        reference: reference transcripts fasta file.
        output: output bam file.
    """
    tool = "minimap2"
    parms = "--secondary=no --cs -a"
    threads = 4
    minimap2(input, reference, output, tool, parms, threads)
    

def salmon_quantify(input, reference, output):
    """
    Quantify with salmon.
    Args:
        input: input bam file.
        reference: reference transcripts fasta file.
        output: output quantification file.
    """
    
    command = [
        "salmon", "quant", "--noErrorModel", "-p", "4",
        "-t", reference, "-l", "SF", 
        "-a", input, "-o", output
        ]
    run_command(command)

def matrix_generate(input_dirs, output, count_type):
    """
    Generate count matrix from salmon results.
    Args:
        input_dirs:  Regular matching pattern for salmon results directories, such as 'path/to/*_quant'.
        output: output file name.
        count_type: count type to extract ("NumReads" or "TPM").
    """

    if count_type not in ["NumReads", "TPM"]:
        raise ValueError("count_type should be 'NumReads' or 'TPM'")
    
    count_data = OrderedDict()
    transcript_ids = None

    if isinstance(input_dirs, str):
        input_dirs = sorted(glob.glob(input_dirs))
    elif not input_dirs:
        raise ValueError("input_dirs should not be empty")
    
    for sample_dir in input_dirs:
        sample_name = os.path.basename(os.path.normpath(sample_dir))
        quant_file = os.path.join(sample_dir, "quant.sf")
        
        if not os.path.exists(quant_file):
            raise FileNotFoundError(f"quant.sf not found in {sample_dir}")
        
        df = pd.read_csv(quant_file, sep='\t')
        
        if "Name" not in df.columns or count_type not in df.columns:
            raise ValueError(f"Required columns missing in {quant_file}")
        
        if transcript_ids is None:
            transcript_ids = df["Name"].tolist()
            count_data["transcript_id"] = transcript_ids
        else:
            if transcript_ids != df["Name"].tolist():
                raise ValueError("Transcript IDs mismatch between samples. "
                                 "Ensure all samples used the same reference.")
        count_data[sample_name] = df[count_type].tolist()
    
    count_df = pd.DataFrame(count_data)
    count_df.set_index("transcript_id", inplace=True)
    count_df.to_csv(output, sep='\t')