#!/usr/bin/env bash
cd /gpfs/data01/bennerlab/home/sjroth/ReadThrough/Benchmarking
for i in {1..10}
do
    ARTDeco -bam-files-dir aligned_files/ -gtf-file modified_genes.gtf \
    -chrom-sizes-file ../../genomes/hg38.chrom.sizes -cpu 50 -layout SE -stranded True -orientation Reverse \
    -skip-bam-summary -meta-file meta.txt -comparisons-file comparisons.txt
    rm -r diff_exp* dogs/ preprocess_files/ quantification/ readthrough/ summary_files/
done