# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250217
Description:
function map.
"""
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
        dupfile = {}
        for line in file:
            count = 1
            fields = line.strip().split("\t")
            if dupfile.get(fields[4]) is None:
                dupfile[fields[4]] = fields + [str(count)]
            else:
                count += int(dupfile[fields[4]][-1])
                dupfile[fields[4]][-1] = str(count)
    
        for key, values in dupfile.items():
            output_filename = prefix + "_" + values[6] + ".bed" 
            with open(output_filename, 'a') as output_file:
                output_file.write("\t".join(values) + "\n")
