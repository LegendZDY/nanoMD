import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from ..utils.modifications import form_reads_get_modifications

app = typer.Typer()

@app.command()
def detectMod(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    sam: Annotated[str, typer.Option("--sam", "-s", help="mapping sam/bam file.")],
    bed: Annotated[str, typer.Option("--bed", "-b", help="bed file for modification sites.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output file.")],
    pvalue: Annotated[float, typer.Option("--pvalue", "-p", help="pvalue cutoff for modification sites.")]=0.98,
    ):
    """
    detect modification sites in input fastq files.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="detecting modification sites...", total=None)
            start=time.time()
            mod = form_reads_get_modifications(input, sam, bed, output, pvalue)
            mod.get_mod_position_with_sam()
            end=time.time()
            time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"detecting modification Done, time cost: {time_cost}")
            progress.add_task(description=f"detecting modification sites Done, time cost: {time_cost}", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="detecting modification sites Failed", total=None)
            exit(1)