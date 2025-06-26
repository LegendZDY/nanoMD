import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import check_path_exists
from ..utils.quantify import matrix_generate

app = typer.Typer()

@app.command()
def matrix(
    input_dirs: Annotated[str, typer.Option("--input_dirs", "-i", help="Regular matching pattern for salmon results directories, such as 'path/to/*_quant'.")],
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

            output_count = "matrix_count.tsv"
            if not check_path_exists(output_count):
                matrix_generate(input_dirs, output=output_count, count_type="NumReads")
            output_tpm = "matrix_tpm.tsv"
            if not check_path_exists(output_tpm):
                matrix_generate(input_dirs, output=output_tpm, count_type="TPM")

            end=time.time()
            time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"Generate count matrix Done, time cost: {time_cost}")
            progress.add_task(description=f"Generate count matrix Done, time cost: {time_cost}", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="Generate count matrix Failed", total=None)
            exit(1)