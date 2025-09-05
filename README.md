# README

This project is an effort to learn python by solving [bioinformatic challanges](https://rosalind.info/problems/locations/).

## Environment setup

> [!NOTE]
> Assumes the packages manager Conda is already installed on the system

Python is ran inside a conda environment with all the nessessary dependencies installed and contained within. The environment can be setup in multiple ways, here we build the environment from a single file. Paste the following into `environment.yml`

```yml
name: bioinformatics
channels:
  - conda-forge
dependencies:
  - python=3.13
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


---


+ env: conda dependencies
+ webscrape code explaination
+ statistics for prediction





Try the script. I am using marimo as notebook.

```sh
marimo edit stryket.py
```