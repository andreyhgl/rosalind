# README

This project is an effort to learn python by solving [bioinformatic challenges](https://rosalind.info/problems/locations/).

<details><summary>Environment setup</summary>

  ## Environment setup
  
  > [!NOTE]
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

```sh
python bin/ini.py data/rosalind_ini.txt
```

