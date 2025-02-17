# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250217
Description:
function map.
"""
import subprocess

def split_mod(input, prefix):
    """
    split modification bed file into single bed files.
    Args:
        input: input modification bed file.
        prefix: output file prefix.
    Returns:
        None.
    """
    with open(input, 'r') as file:
        for line in file:
            fields = line.strip().split("\t")
            if len(fields) >= 7:
                output_filename = prefix + "_" + fields[6] + ".bed"
                with open(output_filename, 'a') as output_file:
                    output_file.write(line)