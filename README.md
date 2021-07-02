# ELeFHAntdev
Ensemble Learning for Harmonization and Annotation of Single Cells (ELeFHAnt) provides an easy to use R package for users to annotate clusters of single cells and harmonize labels across single cell datasets to generate a unified atlas. It provides users with a flexibility of choosing a machine learning based classifiers or let ELeFHAnt automatically use the power of robust classifiers like randomForest and SVM (Support Vector Machines) to make predictions. It has two functions 1) CelltypeAnnotation 2) LabelHarmonization.

## Active development git repository for ELeFHAnt
Code is under active development. This is can be used for test purposes or as beta-version. For stable release visit: https://github.com/praneet1988/ELeFHAnt

# Installation
```
library(devtools)
devtools::install_github('praneet1988/ELeFHAntdev')
library(ELeFHAntdev)
```
