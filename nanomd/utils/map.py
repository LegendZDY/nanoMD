# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250114
Description:
function map.
"""
import subprocess

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
    elif suffix == "bam":
        command = [
            tool, *params_list, "-t", str(threads),
            reference, input, "|", "samtools", 
            "sort", "-@", str(threads), "-o", output, 
            "&&", "samtools", "index", output
            ]
    try:
        print(f"Running command: {' '.join(command)}")
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        exit(1)