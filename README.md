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

## NCBI's GenBank query

The NCBI GenBank contains all annotated DNA sequences, with their transcripts and proteins. To extract entries from this database, use NCBI search engine [Entrez](https://www.ncbi.nlm.nih.gov/search/). [Biopython](https://biopython.org/) is a python library with biological computational tools.

<details>
  <summary>Code</summary>

  ```py
  #!/usr/bin/env python
  
  from Bio import Entrez
  ```

</details>