# README

This project is an effort to learn python by solving [bioinformatic challenges](https://rosalind.info/problems/locations/).

<details>
  <summary>Environment setup</summary><br>
  
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

To read a file, use the function `open`.<br>
Add the statement `with` to close the file after read.<br>
Use the `read()` method for `open` to read the content.

```py
def read_file(file):
  with open(file, "r") as f:

    # .strip() drops the last white space
    content = f.read().strip()
  return content
```

Count the nucleotides in a DNA string. Instead of hard-coding the nucleotides lets pick out any unique character w/ the function `set()`.

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

Put it all together

```sh
python bin/ini.py data/rosalind_ini.txt
```