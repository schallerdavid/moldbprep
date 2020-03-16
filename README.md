moldbprep
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/schallerdavid/moldbprep.svg?branch=master)](https://travis-ci.com/schallerdavid/moldbprep)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/schallerdavid/moldbprep/branch/master?svg=true)](https://ci.appveyor.com/project/schallerdavid/moldbprep/branch/master)
[![codecov](https://codecov.io/gh/schallerdavid/moldbprep/branch/master/graph/badge.svg)](https://codecov.io/gh/schallerdavid/moldbprep/branch/master)

Prepare, standardize and merge molecule databases for virtual screening.

## Install

#### Clone this repository

Open a new terminal and clone this repository
```bash
cd ~
git clone https://github.com/schallerdavid/moldbprep
```

#### Install dependencies

PyRod is written in python 3.8 and uses molvs (>=0.1.1, is shipped with RDKit and Pandas), RDKit and Pandas which can be easily installed using conda:

```bash
conda create -n moldbprep -c conda-forge molvs python=3.8
```

#### Create alias for your bash

```bash
echo 'alias moldbprep="python3 ~/moldbprep/moldbprep.py"' >> ~/.bashrc
```

## Run examples

#### Load conda environment

Activate conda environment.
```bash
source activate moldbprep
```

Run standardization of sample dbs.
```bash
moldbprep -i /home/david/moldbprep/moldbprep/data/db1.sdf,/home/david/moldbprep/moldbprep/data/db2.sdf,/home/david/moldbprep/moldbprep/data/db3.sdf -o /home/david/Documents/moldbprep -p 4
```

## Run tests

Activate conda environment.
```bash
source activate moldbprep
```

Install pytest
```bash
conda install pytest
```

Run tests

```bash
cd moldbprep
pytest
```

## Copyright

Copyright (c) 2020, David Schaller


### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
