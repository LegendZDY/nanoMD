import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from ..utils.map import minimap2map

app = typer.Typer()

@app.command()
def isoform(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    reference: Annotated[str, typer.Option("--reference", "-r", help="reference transcripts path.")],
    output: Annotated[str, typer.Option("--output", "-o", help="output for output sam/bam files.")],
    tool: Annotated[str, typer.Option("--tool", help="minimap2.")]="minimap2",
    parms: Annotated[str, typer.Option("--parms", help="minimap2 parameters for mapping.")]="--secondary=no --cs -a",
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads.")]=4,
    ):
    """
    Mapping of nanopore reads to transcripts reference.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="map transcripts reference...", total=None)
            start=time.time()
            minimap2map(input, reference, output, tool, parms, threads)
            end=time.time()
            time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"map transcripts reference Done, time cost: {time_cost}")
            progress.add_task(description=f"map transcripts reference Done, time cost: {time_cost}", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="map transcripts reference Failed", total=None)
            exit(1)
    