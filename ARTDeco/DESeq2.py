'''
Module that contains functions for running basic DESeq2 operations. Used for both diff_exp_read_in and diff_exp_dogs
modes.
'''
import pandas as pd
import os
from itertools import combinations
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri,Formula
pandas2ri.activate()
from rpy2.robjects.packages import importr
deseq = importr('DESeq2')

'''
Define a function that can reformat a meta file.
'''
def reformat_meta(meta_file,out_dir):

    meta_df = pd.read_csv(meta_file,sep='\t')
    meta_df.Experiment = meta_df.Experiment.str.replace('-','_')
    meta_df.Experiment = meta_df.Experiment.str.replace(' ','_')
    meta_df.Group = meta_df.Group.str.replace('-','_')
    meta_df.Group = meta_df.Group.str.replace(' ','_')
    meta_df = meta_df.sort_values('Experiment')
    meta_df.to_csv(os.path.join(out_dir,'meta.reformatted.txt'),sep='\t',index=False)

'''
Define a function that can reformat a comparisons file.
'''
def reformat_comparisons(comparisons_file,out_dir):

    comparisons_df = pd.read_csv(comparisons_file,sep='\t',header=None)
    comparisons_df[0] = comparisons_df[0].str.replace('-','_')
    comparisons_df[0] = comparisons_df[0].str.replace(' ','_')
    comparisons_df[1] = comparisons_df[1].str.replace('-','_')
    comparisons_df[1] = comparisons_df[1].str.replace(' ','_')
    comparisons_df.to_csv(os.path.join(out_dir,'comparisons.reformatted.txt'),sep='\t',index=False,header=False)

'''
Define a function that can generate an all-by-all comparisons file.
'''
def generate_comparisons(meta_file,out_dir):

    meta_df = pd.read_csv(meta_file,sep='\t')
    f = open(os.path.join(out_dir,'comparisons.reformatted.txt'),'w')
    for combo in combinations(meta_df.Group.unique(),2):
        f.write('\t'.join(combo)+'\n')
    f.close()

'''
Define a function that can load a DESeq2 dataset.
'''
def load_deseq_dataset(count_file,meta_file):

    #Load count dataframe.
    count_df = pd.read_csv(count_file,sep='\t',index_col=0)
    del count_df['Length']
    count_df = count_df.astype(int)
    count_r_df = pandas2ri.py2ri(count_df)

    #Load meta dataframe.
    meta_r_df = robjects.DataFrame.from_csvfile(meta_file,sep="\t",row_names=1)

    #Set design and load dataset.
    design = Formula('~ Group')
    dds = deseq.DESeqDataSetFromMatrix(countData=count_r_df,colData=meta_r_df,design=design)

    return dds

'''
Define a function that can run DESeq2.
'''
def run_deseq(dds):

    filter_exp = robjects.r('function(x) x[rowSums(counts(x))>10,]')
    dds = filter_exp(dds)
    dds_results = deseq.DESeq(dds)

    return dds_results

'''
Define a function that can return results from a DESeq2 run.
'''
def deseq_results(dds,condition1,condition2,out_dir):

    #Get DESeq2 results.
    to_dataframe = robjects.r('function(x) data.frame(x)')
    res = to_dataframe(deseq.results(dds,contrast=robjects.StrVector(['Group',condition1,condition2])))
    gene_ids = res.rownames
    res = pandas2ri.ri2py_dataframe(res)
    res.index = gene_ids

    #Output results.
    res.to_csv(os.path.join(out_dir,f'{condition1}-{condition2}-results.txt'),sep='\t')