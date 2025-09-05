import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium")


app._unparsable_cell(
    r"""
    https://rosalind.info/problems/locations/
    """,
    name="_"
)


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
    ### [Count unique characters](https://rosalind.info/problems/dna/)

    > The character represents nucleotides in a DNA sequence

    To read a file, use the function `open`.<br>
    Add the statement `with` to close the file after read.<br>
    Use the `read()` method for `open` to read the content.
    """
    )
    return


@app.cell
def _():
    with open("../data/dna.txt", "r") as f:
        dna = f.read().strip()

    print(dna)
    print(dna.count("A"), dna.count("C"), dna.count("G"), dna.count("T"))
    return


if __name__ == "__main__":
    app.run()
