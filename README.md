
# Mouse-GEM: The generic genome-scale metabolic model of _Mus musculus_

[![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2FMouse-GEM.svg)](https://badge.fury.io/gh/sysbiochalmers%2FMouse-GEM)


#### Brief Model Description

This repository contains the latest version of Mouse-GEM, a mouse genome-scale metabolic model.


#### Citation

 > H. Wang, J. L. Robinson, P. Kocabaş, J. Gustafsson, M. Anton, P.-E. Cholley, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. _PNAS_ 118, e2102344118 (2021). [doi.org/10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118)

#### Model Keywords

**Utilisation:** multi-omics integrative analysis, predictive simulation
**Field:** metabolic-network reconstruction  
**Type of Model:** reconstruction, curated  
**Model source:** [Human-GEM](https://doi.org/10.1126/scisignal.aaz1482)   
**Omic source:** genomics; metabolomics   
**Taxonomic name:** _Mus musculus_  
**Taxonomy ID:** [10090](https://identifiers.org/taxonomy:10090)  
**Genome ID:** [GCA_000001635.8](https://identifiers.org/insdc.gca:GCA_000001635.8)  
**Metabolic System:** general metabolism  
**Tissue:**  
**Bioreactor:**    
**Cell type:**  
**Cell line:**  
**Condition:** generic metabolism


### Model Overview

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
|_Mus musculus_ |   Human-GEM |    |  |  |


## Installation

### Required Software
* A functional MATLAB installation (MATLAB 7.3 and higher).
* The [RAVEN toolbox](https://github.com/SysBioChalmers/RAVEN).
* The [COBRA toolbox](https://github.com/opencobra/cobratoolbox) (not necessary for most functionality).


### Dependencies - Recommended Software
* The libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.


### Installation Instructions
* Clone the [main branch](https://github.com/SysBioChalmers/Mouse-GEM/tree/main) of this repository, or [download the latest release](https://github.com/SysBioChalmers/Mouse-GEM/releases/latest).
* Add the directory to your MATLAB path (instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)).


## Usage

#### Loading/saving the model

`Mouse-GEM.mat` (Recommended if on `main` branch)
* Load and save using the built-in MATLAB `load()` and `save()` functions.

`Mouse-GEM.xml` (SBML format)
* Load using the `importModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN))
* Save using the `exportModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN))


## Websites

- [Metabolic Atlas](https://metabolicatlas.org/) enables visualization and exploration of Mouse-GEM content.


### Contributing

Contributions are always welcome! Please read the [contributing guidelines](.github/CONTRIBUTING.md) to get started.
