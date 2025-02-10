import typer
from .modules.mapping import mapping
from .modules.detectMod import detectMod
from .modules.isoformAS import isoformAS

app = typer.Typer(add_completion=False)

@app.callback()
def callback():
    """
    nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A sites, 
    genes, and isoforms, and new mRNA.
    """

app.command(name="mapping")(mapping)
app.command(name="isoformAS")(isoformAS)
app.command(name="detectMod")(detectMod)


if __name__ == "__main__":
    app()