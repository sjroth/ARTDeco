########################################################################################################################
# Script that runs complete DoGFinder pipeline for benchmarking. This script is optimized for multi-threading in order #
# to preserve a fair comparison with ARTDeco. This script should be run in a custom conda environment, which can be    #
# generated using the following commands:                                                                              #
# conda create -n dogfinder python=2.7  openssl=1.0                                                                    #
# conda install samtools=1.3 bedtools=2.25 pysam pybedtools pybigwig rseqc=2.6.4                                       #
# Other installation instructions are in the DoGFinder Github. Uses the same modified GTF file as ARTDeco.             #
########################################################################################################################

import os
import subprocess
import shutil

os.chdir("/gpfs/data01/bennerlab/home/sjroth/ReadThrough/Benchmarking/")

#Load annotation.
os.mkdir('annotationdir')
process = subprocess.Popen(['Get_loci_annotation','-out','annotationdir','-gtf','modified_genes.gtf'])
stdout,stderr = process.communicate()

#Grab bam files.
bam_files = []
for f in os.listdir('aligned_files'):
    if f[-11:] == '.sorted.bam':
        bam_files.append(os.path.join('aligned_files',f))

#Perform DoGFinder's preprocess.
process = subprocess.Popen(['Pre_Process','-Q',str(len(bam_files)),'-bam',','.join(bam_files),'-ref',
                            os.path.join('annotationdir','loci_annotation.bed')])
stdout,stderr = process.communicate()

#Find DoGs. Do this in parallel.
#Grab downsampled BAM files.
ds_bam_files = []
for f in os.listdir('aligned_files'):
    if f[-6:] == 'DS.bam':
        ds_bam_files.append(os.path.join('aligned_files',f))

#Find DoGs in parallel.
os.mkdir('outdir')
processes = []
for f in ds_bam_files:
    process = subprocess.Popen(['Get_DoGs','-out','outdir','-bam',f,'-a',
                                os.path.join('annotationdir','loci_annotation.bed'),'-s','-minDoGLen','4000',
                                '-minDoGCov','0.6','-w','500','-mode','W','-suff',f.split('.')[1]])
    processes.append(process)

for i in range(len(processes)):
    exit_code = processes[i].wait()

#Join all DoGs.
dogs = []
for f in os.listdir('outdir'):
    dogs.append(os.path.join('outdir',f))

process = subprocess.Popen(['Union_DoGs_annotation','-dog',','.join(dogs),'-out','outdir'])
stdout,stderr = process.communicate()

#Get expression for each experiment on its own DoGs.
processes = []
for f in ds_bam_files:
    process = subprocess.Popen(['Get_DoGs_rpkm','-out','outdir','-bam',f,'-s','-dog',
                                os.path.join('outdir','Final_Dog_annotation_'+f.split('.')[1]+'.bed'),'-suff',
                                f.split('.')[1]])
    processes.append(process)

for i in range(len(processes)):
    exit_code = processes[i].wait()

#Get expression for each experiment for the union of DoGs.
processes = []
for f in ds_bam_files:
    process = subprocess.Popen(['Get_DoGs_rpkm','-out','outdir','-bam',f,'-s','-dog',
                                os.path.join('outdir','union_dog_annotation.bed'),'-suff','union_'+f.split('.')[1]])
    processes.append(process)

for i in range(len(processes)):
    exit_code = processes[i].wait()

#Delete all files produced.
shutil.rmtree('annotationdir')
shutil.rmtree('outdir')
for f in os.listdir('aligned_files'):
    if 'DS' in f:
        os.remove(os.path.join('aligned_files',f))