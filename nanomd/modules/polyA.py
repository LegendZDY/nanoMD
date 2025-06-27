import time
import typer
from typing import List
from pathlib import Path
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import check_path_exists
from ..utils.pod5_to_fast5 import convert_to_fast5_with_summary_file

app = typer.Typer()

@app.command()
def polyA(
    inputs: Annotated[Path, typer.Option("--inputs", "-i", help="Dir path of pod5 files.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
    ):
    """
    Detect polyA RNA using mapping sam/bam file and fastq files. And quantify them using salmon.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        progress.add_task(description="Detecting nascent RNA start...", total=None)
        start=time.time()
        
        summary_file=f"{prefix}_summary.csv"
        progress.add_task(description="Getting fatures from sam files...", total=None)
        if not check_path_exists(summary_file):
            convert_to_fast5_with_summary_file([inputs], output, summary_file)
        progress.add_task(description=f"Getting fatures from sam files Done", total=None)

        # output_fastq=f"{output}/{prefix}_nascentRNA.fastq"
        # progress.add_task(description="Getting nascentRNA reads...", total=None)
        # if not check_path_exists(output_fastq):
        #     new_fq(output_file, model, input, output_fastq)
        # progress.add_task(description="Getting nascentRNA reads Done", total=None)

        # output_bam=f"{output}/{prefix}_nascentRNA.bam"
        # progress.add_task(description="Mapping nascentRNA reads to reference...", total=None)
        # if not check_path_exists(output_bam):
        #     salmon_map(output_fastq, reference, output_bam)
        # progress.add_task(description="Mapping nascentRNA reads to reference Done", total=None)

        # output_quant=f"{output}/{prefix}_quant"
        # progress.add_task(description="Quantifying nascentRNA reads ...", total=None)
        # if not check_path_exists(output_quant):
        #     salmon_quantify(output_bam, reference, output_quant)
        # progress.add_task(description="Quantifying nascentRNA reads Done", total=None)

        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detect nascent RNA Done, time cost: {time_cost}")
        progress.add_task(description=f"Detect nascent RNA Done, time cost: {time_cost}", total=None)