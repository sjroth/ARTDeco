########################################################################################################################
# Script that runs the Dogcatcher pipeline with the differential expression module for benchmarking. This script is    #
# optimized for multi-threading in order to preserve a fair comparison with ARTDeco. This script should be run in a    #
# custom conda environment, which can be generated using the following command:                                        #
# conda create -n dogcatcher pandas bioconductor-subread bioconductor-deseq2 bedtools                                  #
# Uses a modified Ensembl GTF. The R template module from the Dogcatcher GitHub has been modified to work with         #
# single-end reads.                                                                                                    #
########################################################################################################################

import os
import subprocess
from multiprocessing import Pool
import shutil

os.chdir("/gpfs/data01/bennerlab/home/sjroth/ReadThrough/Benchmarking/")

software_dir = '/gpfs/data01/bennerlab/home/sjroth/software/Dogcatcher/'

#GTF file.
gtf_file = os.path.join('dogcatcher_gtf','Homo_sapiens.GRCh38.99.withchr.gtf')

#Load annotation.
process = subprocess.Popen(['python',os.path.join(software_dir,'1.0_Dogcatcher_flatten_gtf.py'),
                            '--annotation_file_with_path',gtf_file])
stdout,stderr = process.communicate()

#Discover DoGs.
'''
Define a function that can run DoG discovery for a given experiment.
'''
def discover_dog(args):
    experiment,out_dir = args
    process = subprocess.Popen(['python',os.path.join(software_dir,'2.0_Dogcatcher.py'),'--annotation_file_with_path',
                                gtf_file,'--cpus','25','--BedGraph_input_min_strand',experiment+'_min.bedGraph',
                                '--BedGraph_input_plu_strand',experiment+'_plu.bedGraph','--output_prefix',out_dir,
                                '--window_size','500','--coverage_percentage','60'])
    stdout,stderr = process.communicate()

sample_to_condition = {}
cmds = []
for f in os.listdir('.'):
    if len(f) > 8 and f[-8:] == 'bedGraph':
        sample = f.split('_')[0]
        if sample not in sample_to_condition:
            condition = sample.split('-')[2]+'/'
            sample_to_condition[sample] = condition
            cmds.append((sample,condition))

pool = Pool(processes=2)
pool.map(discover_dog,cmds)
pool.close()

#Find longest contigs for each batch.
condition_to_sample = {}
for sample,condition in sample_to_condition.items():
    if condition in condition_to_sample:
        condition_to_sample[condition].append(sample)
    else:
        condition_to_sample[condition] = [sample]

for condition in condition_to_sample.keys():
    if condition != 'Mock/':
        os.system('cp -r Mock/* '+condition)

mock_experiments = condition_to_sample['Mock/']
processes = []
for condition,samples in condition_to_sample.items():
    if condition != 'Mock/':
        process = subprocess.Popen(['python',os.path.join(software_dir,'2.5_Dogcatcher_filter.py'),'--filter','longest',
                                    '--input_prefix',condition,'--Dogcatcher_plu_strand_list']+
                                   [x+'_plu.bedGraph' for x in mock_experiments+samples]+
                                   ['--Dogcatcher_min_strand_list']+
                                   [x+'_min.bedGraph' for x in mock_experiments+samples]+
                                   ['--output_prefix',condition[:-1]+'-diff/'])
        processes.append(process)

for i in range(len(processes)):
    exit_code = processes[i].wait()

#Run DESeq2.
for condition,samples in condition_to_sample.items():
    if condition != 'Mock/':
        output_dir = condition[:-1]+'-diff/'
        process = subprocess.Popen(['python',os.path.join(software_dir,'3.0_Create_R_subread_DESeq2_script.py'),
                                '--annotation_file_with_path',gtf_file,'--control_BAM_list']+
                                   [os.path.join('aligned_files',x+'.sorted.bam') for x in mock_experiments]+
                                   ['--treatment_BAM_list']+
                                   [os.path.join('aligned_files',x+'.sorted.bam') for x in samples]+
                                   ['--input_R_template_file',os.path.join(software_dir,'R_subread_DEseq2_TEMPLATE.R'),
                                    '--input_prefix',output_dir,'--output_prefix',
                                    os.path.join(output_dir,'initial_Rsubread_DESeq2'),'--cpus','50','--padj','0.05'])
        stdout,stderr = process.communicate()
        process = subprocess.Popen(['R','CMD','BATCH',
                                os.path.join(output_dir,'initial_Rsubread_DESeq2/Rsubread_DESeq2_initial.R'),
                                os.path.join(output_dir,'initial_Rsubread_DESeq2/Rsubread_DESeq2_initial.R.out')])
        stdout,stderr = process.communicate()

#Generate a gtf with read-through and non-significant genes for proper normalization.
for condition in condition_to_sample.keys():
    if condition != 'Mock/':
        output_dir = condition[:-1] + '-diff/'
        process = subprocess.Popen(['python',os.path.join(software_dir,'4.0_Dogcatcher_Rsubread_DESeq2.py'),
                                    '--annotation_file_with_path',gtf_file,'--input_prefix',output_dir,
                                    '--input_prefix_DESeq2',os.path.join(output_dir,'initial_Rsubread_DESeq2'),
                                    '--output_prefix',os.path.join(output_dir,'Dogcatcher_with_non-significant_genes'),
                                    '--padj','0.05'])
        stdout,stderr = process.communicate()
        for dog in ['DOG','ADOG']:
            for strand in ['plu','min']:
                process = subprocess.Popen(['R','CMD','BATCH',
                                    os.path.join(output_dir,'Dogcatcher_with_non-significant_genes',
                                                 f'{strand}_ALL_SAMPLES_{dog}_with_nonsig.R'),
                                    os.path.join(output_dir,'Dogcatcher_with_non-significant_genes',
                                                 f'{strand}_ALL_SAMPLES_{dog}_with_nonsig.R.out')])
                stdout,stderr = process.communicate()

#Filter results.
processes = []
for condition in condition_to_sample.keys():
    if condition != 'Mock/':
        output_dir = condition[:-1]+'-diff/'
        process = subprocess.Popen(['python',os.path.join(software_dir,'5.0_filter_sig_DESeq2.py'),
                                    '--annotation_file_with_path',gtf_file,'--input_prefix',output_dir,
                                    '--output_prefix',os.path.join(output_dir,'FINAL_OUT'),
                                    '--input_DOG_DESeq2_plu_sense_file',
                                    os.path.join(output_dir,'Dogcatcher_with_non-significant_genes',
                                                 'plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv'),
                                    '--input_DOG_DESeq2_min_sense_file',
                                    os.path.join(output_dir,'Dogcatcher_with_non-significant_genes',
                                                 'min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv'),
                                    '--input_ADOG_DESeq2_plu_antisense_file',
                                    os.path.join(output_dir,'Dogcatcher_with_non-significant_genes',
                                                 'plu_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv'),
                                    '--input_ADOG_DESeq2_min_antisense_file',
                                    os.path.join(output_dir,'Dogcatcher_with_non-significant_genes',
                                                 'min_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv'),
                                    '--padj','0.05'])
        processes.append(process)

for i in range(len(processes)):
    exit_code = processes[i].wait()

#Remove all files generated by script.
#Gene annotation files.
for f in os.listdir('dogcatcher_gtf'):
    fname = os.path.join('dogcatcher_gtf',f)
    if fname != gtf_file:
        os.remove(fname)

#DoG files.
for sample in sample_to_condition.keys():
    shutil.rmtree(sample)
for condition in set(sample_to_condition.values()):
    shutil.rmtree(condition)
    if condition != 'Mock/':
        shutil.rmtree(condition[:-1]+'-diff/')