#!/usr/bin/env bash
for i in {1..10}
do
    snakemake -s DoGFinder -j 50
    rm /gpfs/data01/bennerlab/home/sjroth/ReadThrough/Benchmarking/aligned_files/*sorted*
done