import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from ..utils.map import minimap2map

app = typer.Typer()

@app.command()
def gene(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    reference: Annotated[str, typer.Option("--reference", "-r", help="reference genome path.")],
    output: Annotated[str, typer.Option("--output", "-o", help="output for output sam/bam files.")],
    tool: Annotated[str, typer.Option("--tool", help="minimap2.")]="minimap2",
    parms: Annotated[str, typer.Option("--parms", help="minimap2 parameters for mapping.")]="--secondary=no --cs -ax",
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads.")]=4,
    ):
    """
    Mapping of nanopore reads to a reference genome.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="map reference...", total=None)
            start=time.time()
            minimap2map(input, reference, output, tool, parms, threads)
            end=time.time()
            progress.add_task(description=f"map reference Done, time cost: {end-start}s", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="map reference Failed", total=None)
            exit(1)
    