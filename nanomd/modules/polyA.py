import time
import typer, glob
from pathlib import Path
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import check_path_exists, minimap2
from ..utils.polyAtools import convert_to_fast5_with_summary_file, index_fastq, detect_polyA, PolyADetector

app = typer.Typer()

@app.command()
def polyA(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq file.")],
    transcriptome: Annotated[str, typer.Option("--transcriptome", "-f", help="Reference transcriptome fasta file path.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output file path.")],
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")],
    min_a_length: Annotated[int, typer.Option("--min-a-length", "-a", help="Minimum length of polyA tail.")]=6,
    max_non_a: Annotated[int, typer.Option("--max-non-a", "-n", help="Maximum number of non-A characters in polyA tail.")]=3,
    pod5s: Annotated[str, typer.Option("--pod5s", help="Regular matching pattern for pod5 files, such as 'path/to/*pod5'.")]=None,
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads to use.")]=8,
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
        
        if pod5s != None:
            summary_file=f"{prefix}_summary.txt"
            input_name=Path(input).name
            output_fast5=f"{output}/fast5"
            if isinstance(pod5s, str):
                pod5s_dir = sorted([Path(p) for p in glob.glob(pod5s)])
            elif not pod5s:
                raise ValueError("pod5s should not be empty")
            progress.add_task(description="Pod5 to fast5...", total=None)
            if not check_path_exists(summary_file):
                convert_to_fast5_with_summary_file(pod5s_dir, output_fast5, summary_file, input_name)
            progress.add_task(description=f"Pod5 to fast5 Done", total=None)

            output_index=f"{input}.index"
            progress.add_task(description="Indexing reads...", total=None)
            if not check_path_exists(output_index):
                index_fastq(output_fast5, summary_file, input)
            progress.add_task(description="Indexing reads Done", total=None)

            output_ploya=f"{output}/{prefix}_polyA.tsv"
            sort_bam=f"{output}/{prefix}_polyA.sorted.bam"
            progress.add_task(description="Detecting polyA...", total=None)
            progress.add_task(description="Mapping reads to transcriptome...", total=None)
            if not check_path_exists(sort_bam):
                minimap2(input, transcriptome, sort_bam, params="-ax map-ont", threads=threads)
            progress.add_task(description="Mapping reads to transcriptome Done", total=None)
            if not check_path_exists(output_ploya):
                detect_polyA(input, sort_bam, transcriptome, output_ploya, threads=threads)
            progress.add_task(description="Detecting polyA Done", total=None)
        else:
            output_ploya=f"{output}/{prefix}_polyA.tsv"
            sort_bam=f"{output}/{prefix}_polyA.sorted.bam"
            progress.add_task(description="Detecting polyA...", total=None)
            progress.add_task(description="Mapping reads to transcriptome...", total=None)
            if not check_path_exists(sort_bam):
                minimap2(input, transcriptome, sort_bam, params="-ax map-ont", threads=threads)
            progress.add_task(description="Mapping reads to transcriptome Done", total=None)
            if not check_path_exists(output_ploya):
                polya = PolyADetector(sort_bam, output_ploya, min_a_length, max_non_a)
                polya.analyze()
            progress.add_task(description="Detecting polyA Done", total=None)


        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detect polyA Done, time cost: {time_cost}")
        progress.add_task(description=f"Detect polyA Done, time cost: {time_cost}", total=None)