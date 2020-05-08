'''
Module that contains functions for intergenic mode.
'''
import subprocess
import os
from multiprocessing import Pool
from .misc import load_exp
import functools
import pandas as pd
import numpy as np

'''
Define a function that can get the gene expression given a tag directory, a GTF file, a normalization method, and 
strandedness.
'''
def get_gene_exp(args):

    tag_directory,gtf_file,norm,stranded,out_file = args

    if stranded:
        strand = ['+']
    else:
        strand = ['both']

    f = open(out_file,'w')
    subprocess.call(['analyzeRepeats.pl',gtf_file,'none',norm,'-strand']+strand+['-count','genes','-d',tag_directory],
                    stdout=f,stderr=subprocess.PIPE)
    f.close()

'''
Define a function that can get the gene expression (both normalized and un-normalized) given a list of tag directories, 
a GTF file, and strandedness. 
'''
def get_multi_gene_exp(tag_dirs,gtf_file,stranded,out_dir,cpu):

    #Format and run commands for getting initial gene expression.
    cmds = []
    for norm in ['-raw','-fpkm']:
        cmds += [(tag_dir,gtf_file,norm,stranded,os.path.join(out_dir,tag_dir.split('/')[-1]+f'.{norm[1:]}.txt'))
                for tag_dir in tag_dirs]

    pool = Pool(processes=min(len(cmds),cpu))
    pool.map(get_gene_exp,cmds)
    pool.close()

    #Join all of these files together.
    raw_dfs = []
    fpkm_dfs = []
    for y in [x[-1] for x in cmds]:
        if y[-7:] == 'raw.txt':
            raw_dfs.append(load_exp(y))
            os.remove(y)
        else:
            fpkm_dfs.append(load_exp(y))
            os.remove(y)

    raw_df = functools.reduce(lambda x,y: pd.merge(x,y,on=['ID','Length']),raw_dfs)
    raw_df = raw_df[['ID','Length']+sorted(raw_df.columns[2:])]
    fpkm_df = functools.reduce(lambda x,y: pd.merge(x,y,on=['ID','Length']),fpkm_dfs)
    fpkm_df = fpkm_df[['ID','Length']+sorted(fpkm_df.columns[2:])]

    raw_df.to_csv(os.path.join(out_dir,'gene.exp.raw.txt'),sep='\t',index=False)
    fpkm_df.to_csv(os.path.join(out_dir,'gene.exp.fpkm.txt'),sep='\t',index=False)

'''
Define a function that can get the maximum isoform for all genes when given a gene expression file and a 
gene-to-transcript mapping.
'''
def get_max_isoform(gene_exp_file,gene_to_transcript_file,out_dir):

    #Load gene expression file into dataframe.
    gene_exp = pd.read_csv(gene_exp_file,sep='\t')
    del gene_exp['Length']
    gene_exp = gene_exp.set_index('ID')

    #Get max expression.
    gene_exp['Max Exp'] = gene_exp.max(axis=1)

    #Load gene-to-transcript mapping.
    gene_to_transcript = pd.read_csv(gene_to_transcript_file,sep='\t')

    #Get maximum expression for each gene.
    gene_exp = pd.merge(gene_to_transcript,gene_exp,left_on='Transcript ID',right_index=True)
    max_exp = pd.DataFrame(gene_exp['Max Exp'].groupby(gene_exp['Gene ID']).max())
    max_exp = max_exp.reset_index()

    #Find and store max isoform.
    max_isoform_df = pd.merge(gene_exp,max_exp,left_on=['Gene ID','Max Exp'],right_on=['Gene ID','Max Exp'])
    max_isoform_df = pd.DataFrame(max_isoform_df.groupby('Gene ID').max()['Transcript ID'])
    max_isoform_df = pd.merge(max_isoform_df,gene_to_transcript,on=['Gene ID','Transcript ID'])

    #Output dataframe.
    max_isoform_df.to_csv(os.path.join(out_dir,'max_isoform.txt'),sep='\t',index=False)

'''
Define a function that can get the gene vs. intergenic information as formatted for output.
'''
def get_gene_v_intergenic(gene_count_file,gene_fpkm_file,max_isoform_file,intergenic_count_file,mode,out_file):

    #Load gene expression information for maximum isoform.
    gene_count = pd.read_csv(gene_count_file,sep='\t')
    expts = sorted(gene_count.columns[2:])
    gene_fpkm = pd.read_csv(gene_fpkm_file,sep='\t')
    gene_exp = pd.merge(gene_count,gene_fpkm,on=['ID','Length'],suffixes=(' Gene Count',' Gene FPKM'))
    max_isoform = pd.read_csv(max_isoform_file,sep='\t')
    gene_exp = pd.merge(max_isoform,gene_exp,left_on='Transcript ID',right_on='ID')
    del gene_exp['ID']

    #Load downstream counts.
    intergenic_count = pd.read_csv(intergenic_count_file,sep='\t')

    #Normalize counts by length.
    median_length = np.median(gene_exp.Length.append(intergenic_count.Length))

    for col in gene_exp.columns:
        if 'Count' in col:
            gene_exp[col] = (gene_exp[col]/gene_exp.Length)*median_length
            gene_exp[col[:-10]+'log2 Gene Count'] = np.log2(gene_exp[col]+1)
    del gene_exp['Length']

    for col in intergenic_count.columns[2:]:
        intergenic_count[f'{col} {mode} Count'] = (intergenic_count[col]/intergenic_count.Length)*median_length
        intergenic_count[f'{col} log2 {mode} Count'] = np.log2(intergenic_count[f'{col} {mode} Count']+1)
        del intergenic_count[col]
    del intergenic_count['Length']

    #Join gene and intergenic information.
    gene_w_intergenic = pd.merge(gene_exp,intergenic_count,left_on='Gene ID',right_on='ID')
    del intergenic_count['ID']

    #Get downstream vs. gene log2 ratio.
    for expt in expts:
        gene_w_intergenic[f'{expt} log2Ratio {mode} vs. Gene'] = gene_w_intergenic[f'{expt} log2 {mode} Count']-\
                                                                 gene_w_intergenic[expt+' log2 Gene Count']
        del gene_w_intergenic[f'{expt} log2 {mode} Count']
        del gene_w_intergenic[expt+' log2 Gene Count']

    #Reorder columns for output.
    new_cols = ['Gene ID','Transcript ID']
    for expt in expts:
        new_cols += [expt+' Gene Count',expt+' Gene FPKM',f'{expt} {mode} Count',f'{expt} log2Ratio {mode} vs. Gene']
    gene_w_intergenic = gene_w_intergenic[new_cols]

    gene_w_intergenic.to_csv(out_file,sep='\t',index=False)

'''
Define a function that can deconvolute gene expression using read-in information.
'''
def deconvolute_exp(read_in_file,out_file):

    #Load read-in information.
    read_in = pd.read_csv(read_in_file, sep='\t')

    #Get corrected expression.
    expts = sorted(set([col.split()[0] for col in read_in.columns[2:]]))
    keep_cols = ['Gene ID', 'Transcript ID']
    for expt in expts:
        read_in[f'{expt} Corrected Count'] = read_in[f'{expt} Gene Count']-read_in[f'{expt} Read-In Count']
        read_in[f'{expt} Corrected Count'].clip(lower=0,inplace=True)
        keep_cols += [f'{expt} Gene Count',f'{expt} Corrected Count']
    read_in = read_in[keep_cols]

    #Output corrected expression.
    read_in.to_csv(out_file,sep='\t',index=False)

'''
Define a function that can take in a read-in vs. gene expression file and output a dataframe of primary induction and
read-in gene assignments.
'''
def assign_genes(read_in_file,read_in_threshold,read_in_fpkm,out_file,gene_types_file,gene_types):

    #Load read-in information.
    read_in_df = pd.read_csv(read_in_file,sep='\t')

    #Load gene type information and limit genes for consideration to those gene types.
    if sum(1 for line in open(gene_types_file)) > 1:

        gene_types_df = pd.read_csv(gene_types_file,sep='\t')

        if gene_types and len(gene_types_df) > 0:
            gene_types_df = gene_types_df[gene_types_df['Gene Type'].isin(gene_types)]
            read_in_df = read_in_df[read_in_df['Gene ID'].isin(gene_types_df['Gene ID'])]

    #Get experiments.
    expts = list(sorted(set([x.split()[0] for x in read_in_df.columns[2:]])))

    #For each experiment, infer read-in genes.
    expt_dfs = []
    for expt in expts:

        #Get columns for experiment.
        cols = ['Gene ID','Transcript ID',expt+' Gene FPKM',expt+' log2Ratio Read-In vs. Gene']
        expt_df = read_in_df[cols].copy()

        #Assign read-in genes.
        expt_df[expt+' Assignment'] = 'N/A'
        primary_induction = expt_df[(expt_df[expt+' Gene FPKM'] > read_in_fpkm) &
                                    (expt_df[expt+' log2Ratio Read-In vs. Gene'] < read_in_threshold)].index
        read_in = expt_df[(expt_df[expt+' Gene FPKM'] > read_in_fpkm) &
                          (expt_df[expt+' log2Ratio Read-In vs. Gene'] > read_in_threshold)].index
        expt_df.loc[primary_induction,expt+' Assignment'] = 'Activated'
        expt_df.loc[read_in,expt+' Assignment'] = 'Read-In'
        del expt_df[expt+' Gene FPKM']

        expt_dfs.append(expt_df)

    #Join these dataframes.
    inference_df = functools.reduce(lambda x,y: pd.merge(x,y,on=['Gene ID','Transcript ID']),expt_dfs)
    inference_df.to_csv(out_file,sep='\t',index=False)

'''
Define a function that can summarize read-in and readthrough levels.
'''
def summarize_readthrough_stats(readthrough_file,expts,mode,num_genes,gene_types_file,gene_types):

    df = pd.read_csv(readthrough_file,sep='\t')
    if sum(1 for line in open(gene_types_file)) > 1:
        gene_types_df = pd.read_csv(gene_types_file,sep='\t')

        if gene_types and len(gene_types_df) > 0:
            gene_types_df = gene_types_df[gene_types_df['Gene Type'].isin(gene_types)]
            df = df[df['Gene ID'].isin(gene_types_df['Gene ID'])]

    summary_dfs = []
    output = mode+' Summary\n'
    for expt in expts:
        summary = df[[expt+' Gene FPKM',f'{expt} log2Ratio {mode} vs. Gene']]
        summary = summary.nlargest(num_genes,expt+' Gene FPKM')
        summary = summary[f'{expt} log2Ratio {mode} vs. Gene'].describe()
        summary_dfs.append(summary)

    output_df = functools.reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True),summary_dfs)

    return output+'\n'.join(output_df.to_string().split('\n'))

'''
Define a function that can summarize read-in assignments.
'''
def summarize_read_in_assignments(assignment_file,expts,read_in_threshold,read_in_fpkm):

    df = pd.read_csv(assignment_file,sep='\t')

    summary_dfs = []
    output = f'Read-In Assignments for each experiment for threshold of {read_in_threshold} and FPKM >= '+\
             f'{read_in_fpkm}\n'
    for expt in expts:
        summary = pd.DataFrame(df.groupby(expt+' Assignment').count()['Gene ID'])
        summary.columns = [f'{summary.index.name}']
        summary.index.name = 'Assignment'
        summary_dfs.append(summary)

    output_df = functools.reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True),summary_dfs)

    return output+'\n'.join(output_df.to_string().split('\n'))