# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250114
Description:
function map.
"""
import time
import subprocess
from rich.progress import Progress, SpinnerColumn, TextColumn


def minimap2map(input, reference, prefix, tool, params, threads):
    """
    map with minimap2.
    Args:
        input: input fasta file.
        reference: reference fasta file.
        prefix: output prefix.
        threads: number of threads.
        tool: minimap2 tool path.
        params: minimap2 params.
    """
    command = [
        tool, params, 
        "-t", threads,
        reference, input,
        "-o", f"{prefix}.sam"
    ]

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="minimap2 mapping...", total=None)
            start=time.time()
            subprocess.run(command, check=True)
            end=time.time()
            progress.add_task(description=f"minimap2 mapping Done, time cost: {end-start}s", total=None)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
            progress.add_task(description="minimap2 mapping Failed", total=None)
            exit(1)