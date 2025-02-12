import time
import typer
from rich.progress import Progress, SpinnerColumn, TextColumn
from typing_extensions import Annotated
from ..utils.modifications import form_reads_get_modifications

app = typer.Typer()

@app.command()
def detectMod(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    sam: Annotated[str, typer.Option("--sam", "-a", help="mapping sam file.")],
    bed: Annotated[str, typer.Option("--bed", "-b", help="bed file for modification sites.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output file.")],
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
            mod = form_reads_get_modifications(input, sam, bed, output)
            mod.get_mod_position_with_sam()
            end=time.time()
            progress.add_task(description=f"detecting modification sites Done, time cost: {end-start}s", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="detecting modification sites Failed", total=None)
            exit(1)
    
    

