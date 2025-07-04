import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import check_path_exists
from ..utils.quantify import matrix_generate
from ..utils.polyAtools import polya_matrix_generate

app = typer.Typer()

@app.command()
def matrix(
    inputs: Annotated[str, typer.Option("--inputs", "-i", help="Regular matching pattern for salmon results directories or polyA files, such as 'path/to/*_quant' or 'path/to/*polyA.tsv'.")],
    control_names: Annotated[str, typer.Option("--control_names", "-c", help="Control sample dir name or polyA file name separated by comma, such as 'control1_quant,control2_quant,control3_quant' or 'control1_polyA.tsv,control2_polyA.tsv,control3_polyA.tsv'.")],
    type: Annotated[str, typer.Option("--type", "-t", help="Input count file type, either 'salmon' or 'polyA'.")] = "salmon",
    ):
    """
    Generate count matrix with salmon results.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="Generate count matrix...", total=None)
            start=time.time()
            if type == "salmon":
                output_count = "matrix_count.tsv"
                if not check_path_exists(output_count):
                    matrix_generate(input_dirs=inputs, control_dirs_names=control_names, output=output_count,
                                    count_type="NumReads")
                output_tpm = "matrix_tpm.tsv"
                if not check_path_exists(output_tpm):
                    matrix_generate(input_dirs=inputs, control_dirs_names=control_names, output=output_tpm,
                                    count_type="TPM")
            elif type == "polyA":
                output_count = "matrix_polyA.tsv"
                if not check_path_exists(output_count):
                    polya_matrix_generate(input_files=inputs, control_filess_names=control_names, output_file=output_count)
            else:
                print(f"Error: Unknown count type {type}")

            end=time.time()
            time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"Generate count matrix Done, time cost: {time_cost}")
            progress.add_task(description=f"Generate count matrix Done, time cost: {time_cost}", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="Generate count matrix Failed", total=None)
            exit(1)