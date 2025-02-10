import typer
from .modules.gene import gene
from .modules.isoform import isoform
from .modules.detectMod import detectMod
from .modules.isoformAS import isoformAS
from .modules.nascentRNA import nascentRNA

app = typer.Typer(add_completion=False)

@app.callback()
def callback():
    """
    nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A sites, 
    genes, and isoforms, and new mRNA.
    """

app.command(name="gene")(gene)
app.command(name="isoform")(isoform)
app.command(name="isoformAS")(isoformAS)
app.command(name="detectMod")(detectMod)
app.command(name="nascentRNA")(nascentRNA)


if __name__ == "__main__":
    app()