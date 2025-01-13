import typer
import subprocess
from pathlib import Path
from rich.progress import track
from typing_extensions import Annotated
from .utils.modifications import getModifications

app = typer.Typer()

APP_NAME = "nanomd"
app_dir = typer.get_app_dir(APP_NAME)

@app.callback()
def callback():
    """
    nanomd is a package for analyzing nanopore RNA sequencing data.
    """
# 定义一个命令来运行 map.sh 脚本
@app.command()
def detect(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    read1: Annotated[str, typer.Option("--read1", "-1", help="Read1 fastq file.")],
    read2: Annotated[str, typer.Option("--read2", "-2", help="Read2 fastq file.")],
    linker: Annotated[str, typer.Option("--linker", "-l", help="Linker sequence.")],
    rnafq: Annotated[str, typer.Option("--rnafq", "-r", help="Output directory for RNA fastq files.")],
    dnafq: Annotated[str, typer.Option("--dnafq", "-d", help="Output directory for DNA fastq files.")],
    structure: Annotated[str, typer.Option("--structure", "-s", help="RNA-linker-DNA structure.")]="rna-linker-dna",
    threads: Annotated[int, typer.Option(help="Number of threads for QC.")]=4,
    ):
    """
    Split fastq files into RNA and DNA fastq files.
    """
    map_script = Path(__file__).parent / "scripts" / "qc.sh"
    subprocess.run(["bash", str(map_script), input, linker, rnafq, dnafq, structure, str(threads)], check=True)

@app.command()
def map(
    input: Annotated[str, typer.Option(help="Input fastq files.")],
    output: Annotated[str, typer.Option(help="Output directory for mapping results.")],
    reference: Annotated[str, typer.Option(help="Reference genome for mapping.")],
    index: Annotated[str, typer.Option(help="minimap2 index for mapping.")],
    threads: Annotated[int, typer.Option(help="Number of threads for mapping.")]=4,
    ):
    """
    Mapping of nanopore reads to a reference genome.
    """
    map_script = Path(__file__).parent / "scripts" / "map.sh"
    subprocess.run(["bash", str(map_script), input, output, reference, index, str(threads)], check=True)

@app.command()
def pairmerge(
    input: Annotated[str, typer.Option(help="Input sam/bam file.")],
    output: Annotated[str, typer.Option(help="Output directory for pileup results.")],
    threads: Annotated[int, typer.Option(help="Number of threads for pileup.")]=4,
    ):
    """
    Pileup of mapped reads.
    """
    map_script = Path(__file__).parent / "scripts" / "pileup.sh"
    subprocess.run(["bash", str(map_script), input, output, str(threads)], check=True)

@app.command()

def pairmerge(
    pileup: Annotated[str, typer.Option(help="Input pileup file.")],
    rawFa: Annotated[str, typer.Option(help="Input fasta file.")],
    output: Annotated[str, typer.Option(help="Output file for nascentRNA.")],
    hyperTCount: Annotated[int, typer.Option(help="Hyper-T count threshold.")]=1,
    hyperSNP: Annotated[float, typer.Option(help="Hyper-SNP threshold.")]=0.48,
    hyperScoreM: Annotated[int, typer.Option(help="Hyper-ScoreM threshold.")]=-1,
    hyperScoreP: Annotated[int, typer.Option(help="Hyper-ScoreP threshold.")]=1,
    hyperES: Annotated[int, typer.Option(help="Hyper-ES threshold.")]=3,
    hyperNew: Annotated[int, typer.Option(help="Hyper-New threshold.")]=2,
    threads: Annotated[int, typer.Option(help="Number of threads for nascentRNA detection.")]=4,
    ):
    """
    Detection of nascentRNA from mapped reads.
    """
    map_script = Path(__file__).parent / "scripts" / "nascentRNAdetect.py"
    subprocess.run(["python", str(map_script), pileup, rawFa, output, str(hyperTCount), str(hyperSNP), str(hyperScoreM), str(hyperScoreP), str(hyperES), str(hyperNew), str(threads)], check=True)

@app.command()
def rdquantify(
    input: Annotated[str, typer.Option(help="Input sam/bam file.")],
    output: Annotated[str, typer.Option(help="Output directory for quantification results.")],
    group: Annotated[str, typer.Option(help="Group information for quantification.")],
    threads: Annotated[int, typer.Option(help="Number of threads for quantification.")]=4,
    ):
    """
    Quantification of nascentRNA and allRNA.
    """
    map_script = Path(__file__).parent / "scripts" / "quantify.sh"
    subprocess.run(["bash", str(map_script), input, output, group, str(threads)], check=True)
    
@app.command()
def diffexp(
    input: Annotated[str, typer.Option(help="Input directory for quantification results.")],
    output: Annotated[str, typer.Option(help="Output directory for differential expression results.")],
    ):
    """
    Differential expression analysis of nascentRNA and allRNA.
    """
    map_script = Path(__file__).parent / "scripts" / "diffexp.sh"
    subprocess.run(["bash", str(map_script), input, output], check=True)

@app.command()
def altersplice(
    input: Annotated[str, typer.Option(help="Input sam/bam file.")],
    output: Annotated[str, typer.Option(help="Output directory for alternative splicing results.")],
    ):
    """
    Alternative splicing analysis of nascentRNA and allRNA.
    """
    map_script = Path(__file__).parent / "scripts" / "altersplice.sh"
    subprocess.run(["bash", str(map_script), input, output], check=True)

@app.command()
def halflife(
    input: Annotated[str, typer.Option(help="Input directory for quantification results.")],
    output: Annotated[str, typer.Option(help="Output directory for half-life plotting results.")],
    group: Annotated[str, typer.Option(help="Group information for half-life plotting.")],
    type: Annotated[str, typer.Option(help="Type of gene for half-life plotting.")]="gene",
    ):
    """
    Half-life plotting of gene and isoform.
    """
    plot_script = Path(__file__).parent / "scripts" / "halflife.R"
    subprocess.run(["Rscript", str(plot_script), input, output, group, type], check=True)

