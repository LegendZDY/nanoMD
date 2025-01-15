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


def minimap2map(inputfasq, reference, output, threads, soft="minimap2", params="--secondary=no --cs -a"):
    """
    A simple implementation of the map function in Python.
    """
    command = [
        soft, params, 
        "-t", threads,
        reference, inputfasq,
        "-o", f"{output}.sam"
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
    
