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
        run_command(command)

    elif suffix == "bam":
        sam_name = output.replace(".bam", ".sam")
        command = [
            tool, *params_list, "-t", str(threads),
            reference, input,  "-o", sam_name
        ]
        run_command(command)
        run_command(["samtools", "view", "-bS", sam_name, "-o", output])
        run_command(["rm", sam_name])
        run_command(["samtools", "index", output])