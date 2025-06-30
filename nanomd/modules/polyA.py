import time
import typer, glob
from pathlib import Path
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import check_path_exists, run_command
from ..utils.polyAtools import convert_to_fast5_with_summary_file, index_fastq, detect_polyA

app = typer.Typer()

@app.command()
def polyA(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq file.")],
    pod5s: Annotated[str, typer.Option("--pod5s", help="Regular matching pattern for pod5 files, such as 'path/to/*pod5'.")],
    transcriptome: Annotated[str, typer.Option("--transcriptome", help="Transcriptome fasta file path.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
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
        sam=f"{output}/{prefix}_polyA.sam"
        bam=f"{output}/{prefix}_polyA.bam"
        sort_bam=f"{output}/{prefix}_polyA.sorted.bam"
        progress.add_task(description="Detecting polyA...", total=None)
        progress.add_task(description="Mapping reads to transcriptome...", total=None)
        if not check_path_exists(sort_bam):
            run_command(f"minimap2 -a -x map-ont {transcriptome} {input} -o {sam}".split())
            run_command(["samtools", "view", "-bS", sam, "-o", bam])
            run_command(["samtools", "sort", "-o", sort_bam, bam])
            run_command(["samtools", "index", sort_bam])
            run_command(["rm", sam, bam])
        progress.add_task(description="Mapping reads to transcriptome Done", total=None)
        if not check_path_exists(output_ploya):
            detect_polyA(input, sort_bam, transcriptome, output_ploya, threads=threads)
        progress.add_task(description="Detecting polyA Done", total=None)

        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detect polyA Done, time cost: {time_cost}")
        progress.add_task(description=f"Detect polyA Done, time cost: {time_cost}", total=None)