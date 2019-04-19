# Automatic Readthrough Transcription DEteCtiOn: ARTDeco
ARTDeco is a pipeline for analyzing and characterizing transcriptional readthrough as described in Roth et al. (2019, in preparation). Broadly, ARTDeco functions to process a set of BAM files such that transcriptional readthrough can be quantified via a variety of measures including read-in levels, readthrough levels, downstream of gene (DoG) transcript detection, and inference of read-in genes.
## Getting started
### Prerequisites
ARTDeco is a Python package that requires Python 3.6 or higher. Additionally, the following packages are required:
```
bedops>=2.4.35
bx-python>=0.8.2
DESeq2>=1.20
Homer>=4.9
numpy>=1.16.2
pandas>=0.24.2
rpy2>=2.9.4
RSeQC>=3.0.0
samtools=1.9
setuptools
```
The easiest way to install these is to use either Anaconda or Miniconda as a package manager using the [Bioconda](https://bioconda.github.io/) channel and the following command:
```
conda install bedops=2.4.* bioconductor-deseq2=1.20.* bx-python=0.8.* homer=4.9.* numpy=1.16.* pandas=0.24.* pybigwig rpy2=2.9.* rseqc=3.0.* samtools=1.9
```
This is recommended because the dependency management of conda is very good. With this in mind, I've created a simple conda environment for running ARTDeco (located in the Conda folder) if you don't feel like adding these packages to your current environment. Go into the Conda directory and run the following code:
```
conda env create -f environment.yml
```
This creates an environment called ARTDeco.

A conda recipe for ARTDeco is in development and should be out fairly soon with Bioconda. It is important to note that Bioconda's version of Homer is a bit buggy so you will not get the full functionality of Homer (which is not necessary here).
### Installation
Once the prerequisites are installed, go into the same directory as setup.py and run the following code:
```
python setup.py install
```
ARTDeco should be installed. In the future, there will be both pip and conda options for installation.
### Preparing files
Once all prerequisites have been installed, there are a few files that are necessary before starting. At a bare minimum, you need a GTF file for your genome of interest as well as a chromosome sizes file. 

ARTDeco converts these GTF files to BED files for manipulation as well as detection of BAM file information (covered in more depth below) using a utility from BEDOPS (gtf2bed). gtf2bed is formatted for GTFs with UCSC standard format and is temperamental with GENCODE files. Run the following command to check if your GTF file works:
``` 
gtf2bed < genes.gtf
```
If the output is the following, then you have an issue:
``` 
Error: Potentially missing gene or transcript ID from GTF attributes (malformed GTF at line [1]?)
```
The solution to this is actually simple:
``` 
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' genes.gtf > modified_genes.gtf
```
Hopefully, BEDOPS will correct this bug or I will write my own parser if possible.

Additionally, ARTDeco requires a chromosome sizes file. This is easily generated using samtools with the following few lines of code:
``` 
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > genome.chrom.sizes
```

There are two other files that need to be generated if the user wants to utilize differential expression in the analyses (covered below). These are the meta file and the comparisons file. The meta file puts labels to each experiment for DESeq2 and the comparisons file tells the differential expression script which comparisons to outcome. These files are tab-delimited and examples are given in the SampleInputs folder in the repository.

The comparisons file is not strictly necessary and ARTDeco will generate one. However, this is not recommended as ARTDeco will do an all-by-all comparisons and this can be time consuming with larger datasets.

ARTDeco will automatically reformat these files when running any differential expression so that they play well with DESeq2. 

The only remaining requirement in terms of files is that you place all of your BAM files in a directory. ARTDeco will generate a directory structure (discussed in the [Interpreting ARTDeco Outputs](#interpreting-artdeco-outputs) section) around these files.
## Running ARTDeco
In this section, standard usage for ARTDeco will be outlined. The next section will explain how to interpret outputs.
## Interpreting ARTDeco Outputs