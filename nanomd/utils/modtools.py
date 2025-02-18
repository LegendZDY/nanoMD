# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250217
Description:
function map.
"""
from collections import defaultdict

def split_mod(input, prefix):
    """
    split modification bed file into single bed files.
    Args:
        input: input modification bed file.
        prefix: output file prefix.
    Returns:
        None.
    """
    filemod = defaultdict(list)
    
    with open(input, 'r') as file:
        for line in file:
            fields = line.strip().split("\t")
            if filemod.get(fields[4]) is None:
                filemod[fields[4]] = fields.append(1)
            else:
                filemod[fields[4]][-1] += 1
    
    for key, value in filemod.items():
        output_filename = prefix + "_" + value[6] + ".bed"
        with open(output_filename, 'a') as output_file:
            output_file.write("\t".join(value))


