import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from ..utils.basetools import check_path_exists
from ..utils.sam_fatures_catch_all import random_forest_data
from ..utils.nascentRNA_fetch import new_fq

app = typer.Typer()

@app.command()
def nascentRNA(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    sam: Annotated[str, typer.Option("--sam", "-s", help="mapping sam/bam file.")],
    base: Annotated[str, typer.Option("--base", "-b", help="choose one of AUGC.")],
    model: Annotated[str, typer.Option("--model", "-m", help="train model.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
    ):
    """
    Detect nascent RNA using mapping sam/bam file and fastq files.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        progress.add_task(description="Detecting nascent RNA start...", total=None)
        start=time.time()
        
        output_file=f"{output}/{prefix}_{base}.csv"
        progress.add_task(description="Getting fatures from sam files...", total=None)

        if not check_path_exists(output_file):
            random_forest_data(sam, base, output_file)

        progress.add_task(description=f"Getting fatures from sam files Done", total=None)

        output_fastq=f"{output}/{prefix}_nascentRNA.fastq"
        progress.add_task(description="Getting nascentRNA reads...", total=None)

        if not check_path_exists(output_fastq):
            new_fq(output_file, model, input, output_fastq)
        progress.add_task(description="Getting nascentRNA reads Done", total=None)
        
        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detect nascent RNA Done, time cost: {time_cost}")
        progress.add_task(description=f"Detect nascent RNA Done, time cost: {time_cost}", total=None)