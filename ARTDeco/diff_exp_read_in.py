'''
Module that contains functions for running diff_exp_read_in mode.
'''
import pandas as pd
import numpy as np
import os
import functools

'''
Define a function that can get all of the condition information.
'''
def get_conditions(conditions_file,meta_file):

    #Get conditions.
    condition1,condition2 = conditions_file.split('/')[-1].split('-')[:2]

    #Load the meta information.
    meta = pd.read_csv(meta_file,sep='\t')
    meta = meta[meta.Group.isin([condition1,condition2])]

    condition1_expts = list(meta[meta.Group == condition1].Experiment)
    condition2_expts = list(meta[meta.Group == condition2].Experiment)

    return condition1,condition2,condition1_expts,condition2_expts

'''
Define a function that can merge read-in information with differential expression output.
'''
def read_in_diff_exp(read_in_file,meta_file,conditions_file,out_dir):

    #Get conditons.
    condition1,condition2,condition1_expts,condition2_expts = get_conditions(conditions_file,meta_file)

    #Load expression data.
    read_in_exp = pd.read_csv(read_in_file,sep='\t')

    #Get average count/FPKM information for given experiments by condition.
    keep_cols = ['Gene ID','Transcript ID']
    for suffix in ['Gene Count','Gene FPKM','Read-In Count']:
        condition1_cols = [f'{expt} {suffix}' for expt in condition1_expts]
        read_in_exp[f'{condition1} Average {suffix}'] = read_in_exp[condition1_cols].mean(axis=1)
        keep_cols.append(f'{condition1} Average {suffix}')

        condition2_cols = [f'{expt} {suffix}' for expt in condition2_expts]
        read_in_exp[f'{condition2} Average {suffix}'] = read_in_exp[condition2_cols].mean(axis=1)
        keep_cols.append(f'{condition2} Average {suffix}')

    read_in_exp = read_in_exp[keep_cols]

    #Get the log2 ratio for the gene vs. read-in.
    read_in_exp[condition1+' Read-In vs. Gene'] = np.log2(read_in_exp[condition1+' Average Read-In Count']+1)-\
                                                  np.log2(read_in_exp[condition1+' Average Gene Count']+1)
    read_in_exp[condition2+' Read-In vs. Gene'] = np.log2(read_in_exp[condition2+' Average Read-In Count']+1)-\
                                                  np.log2(read_in_exp[condition2+' Average Gene Count']+1)

    #Load DESeq2 output.
    diff_exp = pd.read_csv(conditions_file,sep='\t',index_col=0)

    #Merge expression information and DESeq2 output. Format for output.
    read_in_full_info = pd.merge(read_in_exp,diff_exp,left_on='Transcript ID',right_index=True)
    read_in_full_info = read_in_full_info[['Gene ID','Transcript ID','baseMean','log2FoldChange','lfcSE','stat',
                                           'pvalue','padj',condition1+' Average Gene Count',
                                           condition1+' Average Gene FPKM',condition1+' Average Read-In Count',
                                           condition1+' Read-In vs. Gene',condition2+' Average Gene Count',
                                           condition2+' Average Gene FPKM',condition2+' Average Read-In Count',
                                           condition2+' Read-In vs. Gene']]

    #Output merged information.
    read_in_full_info.to_csv(os.path.join(out_dir,f'{condition1}-{condition2}-read_in.txt'),sep='\t',index=False)

'''
Define a function that can assign read-in genes for upregulated genes given thresholds for log2 fold change, p-value,
FPKM, and a read-in threshold.
'''
def assign_read_in_genes(diff_exp_read_in_file,log2FC,pval,FPKM,read_in_threshold,gene_types_file,gene_types,out_dir):

    #Get conditions.
    condition1,condition2 = diff_exp_read_in_file.split('/')[-1].split('-')[:2]

    #Get upregulated genes.
    diff_exp_read_in = pd.read_csv(diff_exp_read_in_file,sep='\t')
    upreg = diff_exp_read_in[(diff_exp_read_in.log2FoldChange > log2FC) & (diff_exp_read_in.padj < pval) &
                             (diff_exp_read_in[condition1+' Average Gene FPKM'] > FPKM)][
        ['Gene ID','Transcript ID','log2FoldChange','padj',condition1+' Read-In vs. Gene']].copy()

    #Make assignments.
    upreg['Assignment'] = 'Primary Induction'
    read_in_index = upreg[upreg[condition1+' Read-In vs. Gene'] > read_in_threshold].index
    upreg.loc[read_in_index,'Assignment'] = 'Read-In'

    #Limit gene types.
    gene_types_df = pd.read_csv(gene_types_file, sep='\t')

    if gene_types and len(gene_types_df) > 0:
        gene_types_df = gene_types_df[gene_types_df['Gene Type'].isin(gene_types)]
        upreg = upreg[upreg['Gene ID'].isin(gene_types_df['Gene ID'])]

    upreg.to_csv(os.path.join(out_dir,f'{condition1}-{condition2}-read_in_assignment.txt'),sep='\t',index=False)

'''
Define a function that can summarize outputs for read-in assignments using differential expression.
'''
def summarize_diff_exp_read_in_assignments(assignment_files,log2FC,pval,read_in_fpkm,read_in_threshold):

    output = 'Summary for read-in assignments using differential expression\nThresholds are log2 fold change > '+\
             f'{log2FC}, p-value < {pval}, FPKM > {read_in_fpkm}, and read-in level threshold is {read_in_threshold}\n'

    summary_dfs = []
    for assignment_file in assignment_files:
        df = pd.read_csv(assignment_file,sep='\t')
        summary = pd.DataFrame(df.groupby('Assignment').count()['Gene ID'])
        summary.columns = [assignment_file.split('/')[-1][:-23]]
        if 'Read-In' not in summary.index:
            summary.loc['Read-In'] = 0
        summary_dfs.append(summary)

    output_df = functools.reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True),summary_dfs)

    return  output+'\n'.join(output_df.to_string().split('\n'))