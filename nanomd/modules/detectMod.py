import typer
import subprocess
from pathlib import Path
from rich.progress import track
from typing_extensions import Annotated
from ..utils.modifications import getModifications

app = typer.Typer()

@app.command()
def detectMod(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    reference: Annotated[str, typer.Option("--reference", "-r", help="reference genome path.")],
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")],
    tool: Annotated[str, typer.Option("--tool", help="minimap2.")]="minimap2",
    parms: Annotated[str, typer.Option("--parms", help="minimap2 parameters for mapping.")]="--secondary=no --cs -a",
    threads: Annotated[int, typer.Option(help="Number of threads.")]=4,
    ):
    """
    detect modification sites in input fastq files.
    """
    getModifications(input, reference, prefix, tool, parms, threads)

