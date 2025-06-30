import time
import typer, glob
from pathlib import Path
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import check_path_exists
from ..utils.polyAtools import convert_to_fast5_with_summary_file, index_fastq, detect_polyA

app = typer.Typer()

@app.command()
def polyA(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq file.")],
    pod5s: Annotated[str, typer.Option("--pod5s", help="Regular matching pattern for pod5 files, such as 'path/to/*pod5'.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
    ):
    """
    Detect polyA with pod5 and fastq files.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        progress.add_task(description="Detect polyA start...", total=None)
        start=time.time()
        
        summary_file=f"{prefix}_summary.txt"
        input_name=Path(input).name
        if isinstance(pod5s, str):
            pod5s_dir = sorted([Path(p) for p in glob.glob(pod5s)])
        elif not pod5s:
            raise ValueError("pod5s should not be empty")
        progress.add_task(description="Pod5 to fast5...", total=None)
        if not check_path_exists(summary_file):
            convert_to_fast5_with_summary_file(pod5s_dir, output, summary_file, input_name)
        progress.add_task(description=f"Pod5 to fast5 Done", total=None)

        # output_fastq=f"{output}/{prefix}_nascentRNA.fastq"
        progress.add_task(description="Indexing reads...", total=None)
        # if not check_path_exists(output_fastq):
        index_fastq(output, summary_file, input)
        progress.add_task(description="Indexing reads Done", total=None)

        # output_bam=f"{output}/{prefix}_nascentRNA.bam"
        # progress.add_task(description="Mapping nascentRNA reads to reference...", total=None)
        # if not check_path_exists(output_bam):
        #     salmon_map(output_fastq, reference, output_bam)
        # progress.add_task(description="Mapping nascentRNA reads to reference Done", total=None)

        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detect polyA Done, time cost: {time_cost}")
        progress.add_task(description=f"Detect polyA Done, time cost: {time_cost}", total=None)