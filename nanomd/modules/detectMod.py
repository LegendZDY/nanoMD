import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from ..utils.modtools import split_mod
from ..utils.basetools import check_path_exists
from ..utils.modifications import form_reads_get_modifications

app = typer.Typer()

@app.command()
def detectMod(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    sam: Annotated[str, typer.Option("--sam", "-s", help="mapping sam/bam file.")],
    bed: Annotated[str, typer.Option("--bed", "-b", help="bed file for transcripts sites.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
    pvalue: Annotated[float, typer.Option("--pvalue", help="pvalue cutoff for modification sites.")]=0.98,
    ):
    """
    detect modification sites in input fastq files.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        progress.add_task(description="Detecting modification start...", total=None)
        start=time.time()
        
        output_file=f"{output}/{prefix}.bed"
        progress.add_task(description="Getting modification from fq files...", total=None)
        if not check_path_exists(output_file):
            mod = form_reads_get_modifications(input, sam, bed, output_file, pvalue)
            mod.get_mod_position_with_sam()
        progress.add_task(description=f"Getting modification from fq files Done", total=None)

        progress.add_task(description="Splitting modification sites...", total=None)
        if not check_path_exists(f"{output}/{prefix}.m6A.bed"):
            split_mod(output_file, prefix)
        progress.add_task(description="Splitting modification sites Done", total=None)

        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detecting modification sites Done, time cost: {time_cost}")
        progress.add_task(description=f"Detecting modification sites Done, time cost: {time_cost}", total=None)