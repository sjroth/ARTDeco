#Automatic Readthrough Transcription DEteCtiOn: ARTDeco
ARTDeco is a pipeline for analyzing and characterizing transcriptional readthrough as described in Roth et al. (2019, in preparation). Broadly, ARTDeco functions to process a set of BAM files such that transcriptional readthrough can be quantified via a variety of measures including read-in levels, readthrough levels, downstream of gene (DoG) transcript detection, and inference of read-in genes.
##Getting started
###Prerequisites
ARTDeco is a Python package that requires Python 3.6 or higher. Additionally, the following packages are required:
```
setuptools
pandas>=0.24.2
rpy2>=2.9.4
numpy>=1.16.2
bx-python>=0.8.2
RSeQC>=3.0.0
bedops>=2.4.35
Homer>=4.9
DESeq2>=1.20
```
The easiest way to install these is to use either Anaconda or Miniconda as a package manager using the [Bioconda](https://bioconda.github.io/) channel and the following command:
```
conda install pandas=0.24.* rpy2=2.9.* numpy=1.16.* bx-python=0.8.* rseqc=3.0.* bedops=2.4.* pybigwig homer=4.9.* bioconductor-deseq2=1.20.*
```
This is recommended because the dependency management of conda is very good. With this in mind, I've created a simple conda environment for running ARTDeco if you don't feel like adding these packages to your current environment. A conda recipe for ARTDeco is in development and should be out fairly soon with Bioconda. It is important to note that Bioconda's version of Homer is a bit buggy so you will not get the full functionality of Homer (which is not necessary here).