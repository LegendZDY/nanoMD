# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250612
Description:
function map.
"""
from basebio import run_command

def minimap2map(input, reference, output, tool, params, threads):
    """
    map with minimap2.
    Args:
        input: input fasta file.
        reference: reference fasta file.
        output: output sam/bam file.
        threads: number of threads.
        tool: minimap2 tool path.
        params: minimap2 params.
    """
    params_list = params.split()
    suffix = output.split(".")[-1]
    if suffix == "sam":
        command = [
            tool, *params_list, "-t", str(threads),
            reference, input, "-o", output
            ]
        use_shell = False
    elif suffix == "bam":
        command = [
            tool, *params_list, "-t", str(threads),
            reference, input, "|", "samtools", 
            "sort", "-@", str(threads), "-o", output, 
            "&&", "samtools", "index", output
            ]
        use_shell = True
    
    run_command(command, use_shell=use_shell)