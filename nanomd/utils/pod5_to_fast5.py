# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250612
Description:
function map.
"""
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List
from more_itertools import chunked

import pod5 as p5
from pod5.tools.utils import DEFAULT_THREADS, collect_inputs, limit_threads
from pod5.tools.pod5_convert_to_fast5 import *

def convert_to_fast5_with_summary_file(
    inputs: List[Path],
    output: Path,
    summary_file: str,
    recursive: bool = False,
    threads: int = DEFAULT_THREADS,
    force_overwrite: bool = False,
    file_read_count: int = 4000,
):
    if output.exists() and not output.is_dir():
        raise FileExistsError("Cannot output to a file")

    threads = limit_threads(threads)

    with ProcessPoolExecutor(max_workers=threads) as executor:
        total_reads = 0
        futures: Dict[Future, Path] = {}

        # Enumerate over input pod5 files
        for input_idx, source in enumerate(
            collect_inputs(inputs, recursive, "*.pod5", threads=threads)
        ):
            # Open the inputs to read the read ids
            with open(summary_file, "w") as fo:
                with p5.Reader(source) as reader:
                    for chunk_idx, read_ids in enumerate(
                        chunked(reader.read_ids, file_read_count)
                    ):
                        dest = (
                            output / f"{source.stem}.{chunk_idx}_{input_idx}.fast5"
                        ).resolve()
                        for read_id in read_ids:
                            fo.write(f"{read_ids}\t{dest}\n")

                        if dest.exists() and not force_overwrite:
                            raise FileExistsError(
                                "Output path points to an existing file and --force-overwrite not set"
                            )

                        kwargs = {
                            "source": source,
                            "dest": dest,
                            "read_ids": read_ids,
                        }
                        futures[executor.submit(convert_pod5_to_fast5, **kwargs)] = dest  # type: ignore

                    total_reads += len(reader.read_ids)

        print(f"Converting pod5s into {len(futures)} fast5 files. Please wait...")

        status = StatusMonitor(file_count=len(inputs))
        status.increment(files_started=len(inputs), read_count=total_reads)

        for idx, future in enumerate(as_completed(futures)):
            (reads_converted, samples_converted) = future.result()

            status.increment(
                files_ended=1,
                sample_count=samples_converted,
                reads_processed=reads_converted,
            )
            status.print_status()

        status.print_status(force=True)

        print("Conversion complete")