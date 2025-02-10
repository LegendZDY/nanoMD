import typer
from .modules.detectmod import app as detectmod
from .modules.mapping import app as mapping

app = typer.Typer(add_completion=False)

@app.callback()
def callback():
    """
    nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A sites, 
    genes, and isoforms, and new mRNA.
    """

app.add_typer(detectmod)
app.add_typer(mapping)

if __name__ == "__main__":
    app()