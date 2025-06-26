# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250612
Description:
function map.
"""
from basebio import run_command
from .map import minimap2map

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
    minimap2map(input, reference, output, tool, parms, threads)
    

def salmon_quantify(input, reference, output):
    """
    map with minimap2.
    Args:
        input: input bam file.
        reference: reference transcripts fasta file.
        output: output quantification file.
    """
    
    command = [
        "salmon", "quant", "-p", "4",
        "-t", reference, "-l", "SF", 
        "-a", input, "-o", output
        ]
    run_command(command)