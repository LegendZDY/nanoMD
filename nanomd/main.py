import typer
from .modules.mapping import app as mapping
from .modules.detectMod import app as detectMod
from .modules.isoformAS import app as isoformAS

app = typer.Typer(add_completion=False)

@app.callback()
def callback():
    """
    nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A sites, 
    genes, and isoforms, and new mRNA.
    """

app.add_typer(mapping, name="mapping", help="Subcommands for isoform analysis.")
app.add_typer(isoformAS, name="isoformAS", help="Subcommands for isoform analysis.")
app.add_typer(detectMod, name="detectMod", help="Subcommands for isoform analysis.")



if __name__ == "__main__":
    app()