# README

This project is an effort to learn python by solving [bioinformatic challenges](https://rosalind.info/problems/locations/).

<details>
  <summary>Environment setup</summary>

  > [!IMPORTANT]
  > Assumes the packages manager Conda is already installed on the system
  
  Python is ran inside a conda environment with all the nessessary dependencies installed and contained within. The environment can be setup in multiple ways, here the environment is built from a single file: `environment.yml`
  
  ```yml
  name: bioinformatics
  channels:
    - conda-forge
  dependencies:
    - python
    - marimo
    - pandas
    - biopython
  ```
  
  Create the environment and "jump" into it
  
  ```sh
  conda env create -f environment.yml -n bioinformatics
  conda activate bioinformatics
  
  # In case new dependancies are needed:
  
  # 1. add them to environmental.yml
  # 2. remove the environment
  #conda env remove -n bioinformatics
  
  # 3. install from file again
  #conda env create -f environment.yml -n bioinformatics
  ```
  
  Start a python notebook (marimo)
  
  ```sh
  marimo edit
  ```
</details>

## [Sequence counting](https://rosalind.info/problems/ini/)

Count number of nucleotides in a given string, read the string from a file.

<details>
  <summary>Code</summary><br>

  To read a file, use the function `open`.<br>
  Add the statement `with` to close the file after read.<br>
  Use the `read()` method for `open` to read the content.<br>
  Wrap in a neat function.
  
  ```py
  def read_file(file):
    with open(file, "r") as f:
  
      # .strip() drops the last white space
      content = f.read().strip()
    return content
  ```
  
  Instead of hard-coding the nucleotides, extract the unique character w/ the function `set()`.<br>
  Use the method `count()` to count the nucleotides.<br>
  Save the number into a string, separate with a _space_.
  
  ```py
  def count_character(content):
    # extract the unique characters from the string, keep in alphabetic order
    chars = "".join(sorted(set(content)))
  
    # assign the counted chars to output
    output = ""
  
    # loop over each char and count, save as string w/ whitespace
    for char in chars:
      output += str(content.count(char)) + " "
    
    print(output.strip())
  ```
  
  Finally, let the script take in an argument for the sequence file, instead of hard-coding the path.
  
  ```py
  import sys
  
  # get first argument
  file = sys.argv[1]
  ```
  
  Put it all together, see [`bin/ini.py`](bin/ini.py)
  
  ```sh
  python bin/ini.py data/rosalind_ini.txt
  ```
</details>

## [NCBI's GenBank query](https://rosalind.info/problems/gbk/)

Given an organism name and two (publication) dates, return the number of counts of nucleotides found in the GenBank database.

> The NCBI GenBank contains all annotated DNA sequences, with their transcripts and proteins. To extract entries from this database, the NCBI search engine [Entrez](https://www.ncbi.nlm.nih.gov/search/) can be used. [Biopython](https://biopython.org/) is a python library with biological computational tools, including search function for Entrez.

<details>
  <summary>Code</summary>

  Maintain two digits for day and month, pad with zero if needed: `2007/2/9` => `2007/02/09`<br>
  Parse the query with the correct quotes

  ```py
  from datetime import datetime
  from Bio import Entrez

  def entrez_search(organism, start_date, end_date):
    # pad dates with zero if needed

    # strptime creates a datetime object
    start_date = datetime.strptime(start_date, "%Y/%m/%d")

    # strftime creates a string
    start_date = start_date.strftime("%Y/%m/%d")

    Entrez.email = "dummy@domain.io"

    # parse the query
    term = '"' + organism + '"' + "[Organism]" + " AND " + '"' + start_date + '"' + "[Publication Date]" + " : " + '"' + end_date + '"' + "[Publication Date]"

    handle = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(handle)
    print(record["Count"])
  ```

  Read the organism and the two dates from file, make variables of them

  ```py
  import sys

  def read_file(file):
    with open(file, "r") as f:
      content = f.read().strip()
    return content

  # get arguments
  file = sys.argv[1]

  # parse file content
  content = read_file(file).split("\n")
  organism, start_date, end_date = content
  ```

  Put it all together, see [`bin/gbk.py`](bin/gbk.py)

  ```sh
  python bin/gbk.py data/rosalind_gbk.txt
  ```
</details>

## 