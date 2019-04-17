'''
Module that contains assorted utils that do not fall under a single usage category. Not including DESeq2 stuff as it
takes time to load.
'''
import subprocess
from multiprocessing import Pool
import pandas as pd
import os

'''
Define a function that can infer experiment format using RSeQC.
'''
def infer_experiment(args):

    #Run infer_experiment.py.
    bam_file,genes_file = args

    p = subprocess.Popen(['infer_experiment.py','-r',genes_file,'-i',bam_file],stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out,err = p.communicate()

    #Parse output.
    out = out.decode("utf-8").strip().split('\n')

    #Infer paired-end or single-end.
    if out[0] == 'This is SingleEnd Data':
        pe = False
    elif out[0] == 'This is PairEnd Data':
        pe = True

    #Infer strandedness.
    stranded = False
    flip = False
    pos = float(out[2].split(':')[1])
    neg = float(out[3].split(':')[1])

    if pos > 0.55 and pos-neg > 0.1:
        stranded = True
    elif neg > 0.55 and pos-neg < -0.1:
        stranded = True
        flip = True

    return bam_file,pe,stranded,flip

'''
Define a function that can infer the experiment type for multiple experiments.
'''
def infer_experiments_group(bam_files,gene_file,cpu):

    cmds = [(bam_file,gene_file) for bam_file in bam_files]
    pool = Pool(processes=cpu)
    formats = pool.map(infer_experiment,cmds)
    pool.close()

    return formats

'''
Define a function that can format a string that describes the format of an inferred output.
'''
def output_inferred_format(inferred_format):

    out_str = ''

    #PE or SE.
    if inferred_format[1]:
        out_str += 'Paired-End'
    else:
        out_str += 'Single-End'

    #Stranded or unstranded.
    if inferred_format[2]:
        out_str += ', Strand-specific, and '
    else:
        out_str += ' and Single-stranded...'

    if inferred_format[3] and inferred_format[2]:
        out_str += 'Reverse-strand oriented...'
    elif not inferred_format[3] and inferred_format[2]:
        out_str += 'Forward-strand oriented...'

    return out_str

'''
Define a function that will load an expression dataframe for a single experiment.
'''
def load_exp(exp_file):

    #Load the file into a dataframe.
    df = pd.read_csv(exp_file,sep='\t')

    #Select relevant columns.
    keep_columns = [df.columns[0]]
    for col in df.columns[1:]:
        if col not in ['Chr','chr','start','end','Strand','strand','Copies','Annotation/Divergence','Peak Score',
                       'Focus Ratio/Region Size']:
            keep_columns.append(col)
    df = df[keep_columns]

    #Rename the ID and tag count columns.
    rename_cols = {df.columns[0]:'ID'}
    tag_cols = []
    for col in df.columns[1:]:
        if col not in ['Start','End','Length']:
            new_name = col.split()[0]
            new_name = new_name.split('/')[-1]
            new_name = new_name.replace('-','_')
            rename_cols[col] = new_name
            tag_cols.append(new_name)
    df = df.rename(rename_cols,axis=1)

    #Create a length column if one does not exist and delete the start and stop columns.
    if 'Length' not in df.columns:
        df['Length'] = df.End-df.Start
        del df['Start']
        del df['End']

    #Reorder columns and return the dataframe.
    tag_cols.sort()
    df = df[['ID','Length']+tag_cols]

    return df

'''
Define a function that can get the expression at intergenic regions given a list of tag directories, a bed file,
strandedness, and normalization.
'''
def get_regions_exp(args):

    tag_directories,region_file,stranded,norm,out_dir,cpu = args

    #Format and run Homer command.
    if stranded:
        strand = ['+']
    else:
        strand = ['both']

    out_file = os.path.join(out_dir,region_file.split('/')[-1][:-3]+f'{norm[1:]}.tmp.txt')
    f = open(out_file,'w')
    subprocess.call(['annotatePeaks.pl',region_file,'none',norm,'-nogene','-noann','-cpu',str(cpu),'-strand']+strand+
                    ['-d']+tag_directories,stdout=f,stderr=subprocess.PIPE)
    f.close()

    #Load into dataframe and output formatted result.
    region_exp_df = load_exp(out_file)
    os.remove(out_file)
    region_exp_df.to_csv(out_file[:-7]+'txt',sep='\t',index=False)