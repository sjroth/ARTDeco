'''
Module that contains functions for get_dogs mode.
'''
from .misc import load_exp,get_regions_exp

import pandas as pd
import os
import subprocess
from multiprocessing import Pool
from bx.intervals.intersection import Intersecter,Interval
import functools

'''
Define a function that can get genes that can be used for DoG discovery.
'''
def get_dog_screening(genes_file,min_len):

    #Load genes.
    genes = pd.read_csv(genes_file,sep='\t',header=None,names=['Chrom','Start','End','Name','Score','Strand'])
    del genes['Score']

    #Separate genes by strand.
    pos_strand = genes[genes.Strand == '+'].copy()
    neg_strand = genes[genes.Strand == '-'].copy()

    #Get test window.
    pos_strand['Downstream Stop'] = pos_strand.End+min_len
    neg_strand['Downstream Stop'] = neg_strand.Start-min_len

    #Limit genes based upon overlaps with minimum length.
    out = []

    #Positive strand.
    pos_strand_vals = pos_strand.values
    out.append(pos_strand_vals[-1][3])

    for i in range(len(pos_strand_vals)-1):

        row = pos_strand_vals[i]
        next_row = pos_strand_vals[i+1]

        #If there's no overlap, add the row to the output.
        if row[0] == next_row[0]:
            if row[5] < next_row[1]:
                out.append(row[3])
        else:
            out.append(row[3])

    #Negative stand.
    neg_strand_vals = neg_strand.values
    out.append(neg_strand_vals[0][3])

    for i in range(1,len(neg_strand_vals)):

        row = neg_strand_vals[i]
        prev_row = neg_strand_vals[i-1]

        #If there's no overlap, add the row to the output.
        if row[0] == prev_row[0]:
            if row[5] > prev_row[2]:
                out.append(row[3])
        else:
            out.append(row[3])

    #Return data frame with only screening genes.
    screening_genes = genes[genes.Name.isin(out)].sort_values(['Chrom','Start'])

    return screening_genes

'''
Define a function that can get intervals for a given gene with a downstream stop and an interval window.
'''
def get_intervals(gene,downstream_stop,window):

    #Generate the intervals.
    out = []
    if gene[4] == '+':
        for interval_start in range(gene[2],downstream_stop-window,window):
            out.append({'Chrom':gene[0],'Start':interval_start,'End':interval_start+window,
                        'Name':f'{gene[3]}-{interval_start}','Score':0,'Strand':gene[4]})
    else:
        for interval_stop in range(gene[1],downstream_stop+window,-window):
            out.append({'Chrom':gene[0],'Start':interval_stop-window,'End':interval_stop,
                        'Name':f'{gene[3]}-{interval_stop}','Score':0,'Strand':gene[4]})

    return out


'''
Define a function that can generate a DoG intervals dataframe and BED file.
'''
def get_all_intervals(args):

    genes,downstream_stop_df,window,out_dir = args

    #Generate all intervals.
    out = []
    downstream_stop_dict = downstream_stop_df.set_index('Name').T.to_dict('list')
    for gene in genes.values:
        out += get_intervals(gene, downstream_stop_dict[gene[3]][0],window)

    #Format and output dataframe.
    intervals_df = pd.DataFrame(out)[['Chrom','Start','End','Name','Score','Strand']].sort_values(['Chrom','Start'])
    intervals_df.to_csv(out_dir,sep='\t',header=False,index=False)


'''
Define a function that can generate the initial screening BED file given a dataframe of screening genes, a minimum length,
and a window size.
'''
def generate_screening_bed(screening_genes,min_len,window,out_dir):

    #Generate initial downstream stop dataframe.
    screening_downstream_stop = []
    for gene in screening_genes.values:
        if gene[4] == '+':
            screening_downstream_stop.append({'Name':gene[3],'Downstream Stop':gene[2]+min_len+window})
        else:
            screening_downstream_stop.append({'Name':gene[3],'Downstream Stop': gene[1]-min_len-window})
    screening_downstream_stop_df = pd.DataFrame(screening_downstream_stop)

    #Generate initial screening dataframe.
    get_all_intervals((screening_genes,screening_downstream_stop_df,window,os.path.join(out_dir,'intervals.bed')))

'''
Define a function that can get the coverage for a DOG screening file given the file, a tag directory and the 
strandedness. Return a list of genes that meet the coverage cutoff.
'''
def get_interval_coverage(args):

    tag_directory,out_dir,bed_file,min_coverage,stranded = args

    #Format and run Homer command.
    if stranded:
        strand = ['+']
    else:
        strand = ['both']

    out_file = os.path.join(out_dir,tag_directory.split('/')[-1]+'.tmp.txt')
    f = open(out_file,'w')
    subprocess.call(['annotatePeaks.pl',bed_file,'none','-fpkm','-nogene','-noann','-strand']+strand+['-d']+
                    [tag_directory],stdout=f,stderr=subprocess.PIPE)
    f.close()

    #Format output.
    df = load_exp(out_file)
    os.remove(out_file)
    df = df[df[df.columns[-1]] > min_coverage]
    del df['Length']
    del df[df.columns[-1]]
    df['Name'] = df.apply(lambda x: x[0].split('-')[0],axis=1)

    return df


'''
Define a function that can get coverage for multiple sets of interval given tag directories, an out directory, a list of
BED files (or a single BED file), the minimum coverage and the strandedness.
'''
def get_multi_interval_coverage(tag_dirs,out_dir,bed_files,min_coverage,stranded,cpu):

    #Screen for coverage threshold.
    cmds = []
    for i in range(len(tag_dirs)):

        if len(bed_files) > 1:
            bed_file = bed_files[i]
        else:
            bed_file = bed_files[0]

        cmds.append((tag_dirs[i],out_dir,bed_file,min_coverage,stranded))

    pool = Pool(processes=min(cpu,len(tag_dirs)))
    screening_coverage_dfs = pool.map(get_interval_coverage,cmds)
    pool.close()

    return screening_coverage_dfs

'''
Define a function that can get a dictionary of downstream genes.
'''
def get_downstream_genes(genes):

    downstream = {}

    #Positive strand case.
    pos_strand_vals = genes[genes.Strand == '+'].values
    for i in range(len(pos_strand_vals)-1):

        gene = pos_strand_vals[i]
        downstream_gene = pos_strand_vals[i+1]

        if gene[0] == downstream_gene[0]:
            downstream[gene[3]] = downstream_gene[3]

    #Negative strand case.
    neg_strand_vals = genes[genes.Strand == '-'].values[::-1]
    for i in range(len(neg_strand_vals)-1):

        gene = neg_strand_vals[i]
        downstream_gene = neg_strand_vals[i+1]

        if gene[0] == downstream_gene[0]:
            downstream[gene[3]] = downstream_gene[3]

    return downstream

'''
Define a function that gets the downstream stop for a given gene when given read-in information in addition to annotation 
and initial screening information.
'''
def get_downstream_stop(gene,genes,downstream,read_in_genes,chrom_sizes):

    #Get downstream gene.
    curr_gene = gene[3]
    if curr_gene in downstream:
        curr_gene = downstream[curr_gene]
        while curr_gene in downstream and curr_gene in read_in_genes:
            curr_gene = downstream[curr_gene]

    #If there is no downstream for the current gene, set the stop to be the end of the chromosome.
    if curr_gene not in downstream:

        #Positive strand case.
        if gene[4] == '+':
            downstream_stop = chrom_sizes[gene[0]]

        #Negative strand case.
        else:
            downstream_stop = 1

    #If the downstream gene is a read-in gene, set the downstream stop to be the stop of that gene.
    elif curr_gene in read_in_genes:

        #Positive strand case.
        if gene[4] == '+':
            downstream_stop = genes[genes.Name == curr_gene].values[0][2]

        #Negative strand case.
        else:
            downstream_stop = genes[genes.Name == curr_gene].values[0][1]

    #Otherwise, set the downstream start to be the start of that gene.
    else:

        #Positive strand case.
        if gene[4] == '+':
            downstream_stop = genes[genes.Name == curr_gene].values[0][1]

        #Negative strand case.
        else:
            downstream_stop = genes[genes.Name == curr_gene].values[0][2]

    return downstream_stop

'''
Define a function that can get the downstream stop for all genes being screened.
'''
def get_all_downstream_stops(args):

    tag_directory,genes_file,screening_genes,read_in_file,chrom_sizes_file = args

    #Get experiment name from tag directory.
    expt_name = tag_directory.split('/')[-1].replace('-','_')

    #Load genes.
    genes = pd.read_csv(genes_file,sep='\t',header=None,names=['Chrom','Start','End','Name','Score','Strand'])
    del genes['Score']

    #Get downstream genes.
    downstream = get_downstream_genes(genes)

    #Load read-in genes.
    read_in_df = pd.read_csv(read_in_file,sep='\t')
    read_in_genes = set(read_in_df[read_in_df[expt_name+' Assignment'] == 'Read-In']['Gene ID'])

    #Load chromosome sizes.
    chrom_sizes = {}
    for line in open(chrom_sizes_file).readlines():
        line = line.strip().split('\t')
        chrom_sizes[line[0]] = int(line[1])

    downstream_stop = []
    for gene in screening_genes.values:
        downstream_stop.append({'Name':gene[3],
                                'Downstream Stop':get_downstream_stop(gene,genes,downstream,read_in_genes,chrom_sizes)})

    return pd.DataFrame(downstream_stop)

'''
Define a function that can generate a full screening bed using a list of tag directories, a genes file, a list of 
screening genes dataframes, a read-in asssignments file, a chromosome sizes file, a window, the number of cpus and an 
output directory.
'''
def generate_full_screening_bed(tag_dirs,genes_file,screening_genes_dfs,read_in_file,chrom_sizes_file,window,cpu,
                                out_dir):

    #Get all of the downstream stop positions based upon read-in information.
    cmds = []
    for tag_dir in tag_dirs:
        cmds.append((tag_dir,genes_file,screening_genes_dfs[tag_dir],read_in_file,chrom_sizes_file))

    pool = Pool(processes=min(cpu,len(tag_dirs)))
    downstream_dfs = pool.map(get_all_downstream_stops,cmds)
    pool.close()

    #Get all of the intervals
    cmds = []
    for i in range(len(tag_dirs)):
        expt_name = tag_dirs[i].split('/')[-1]
        cmds.append((screening_genes_dfs[tag_dirs[i]],downstream_dfs[i],window,os.path.join(out_dir,
                                                                                            f'{expt_name}.bed')))

    pool = Pool(processes=min(cpu,len(tag_dirs)))
    pool.map(get_all_intervals, cmds)
    pool.close()

'''
Define a function that can get DoG coordinates when given coverage information for a given gene, the strand, and the 
window.
'''
def get_dog_coordinates(gene,strand,window):

    gene_vals = gene.values
    start = gene_vals[0][1]
    curr_interval = 1

    if strand == '+':

        while curr_interval < len(gene_vals) and gene_vals[curr_interval][1] == gene_vals[curr_interval-1][1] + window:
            curr_interval += 1

        end = gene_vals[curr_interval-1][1]

        return {'Name':gene_vals[0][0],'Start':start,'End':end+window,'Strand':strand}

    else:

        while curr_interval < len(gene_vals) and gene_vals[curr_interval][1] == gene_vals[curr_interval-1][1]-window:
            curr_interval += 1

        end = gene_vals[curr_interval-1][1]

        return {'Name':gene_vals[0][0],'Start':end-window,'End':start,'Strand':strand}


'''
Define a function that can remove overlaps from a set of DoGs and return the maximum length DoG.
'''
def max_dogs(dogs_df):

    #Generate overlap intervals for each strand.
    pos_strand_intervals = {}
    neg_strand_intervals = {}
    dogs_vals = dogs_df.values
    for dog in dogs_vals:

        if dog[4] == '+':

            if dog[1] not in pos_strand_intervals:
                pos_strand_intervals[dog[1]] = Intersecter()

            pos_strand_intervals[dog[1]].add_interval(Interval(dog[2],dog[3],dog[0]))

        else:

            if dog[1] not in neg_strand_intervals:
                neg_strand_intervals[dog[1]] = Intersecter()

            neg_strand_intervals[dog[1]].add_interval(Interval(dog[2],dog[3],dog[0]))

    #Remove overlaps. Keep DoGs with max coordinates.
    max_dogs = []
    overlapping_genes = set()
    for dog in dogs_vals:

        #Get overlaps.
        if dog[4] == '+':
            overlaps = pos_strand_intervals[dog[1]].find(dog[2],dog[3])
        if dog[4] == '-':
            overlaps = neg_strand_intervals[dog[1]].find(dog[2],dog[3])

        #No overlaps.
        if len(overlaps) == 1:
            name = dog[0]
            start = dog[2]
            end = dog[3]

        #One or more overlaps.
        else:

            #If the gene or an overlap has occurred, skip it.
            if overlaps[0].value in overlapping_genes:
                continue

            #Positive strand case.
            if dog[4] == '+':
                name = overlaps[0].value
                start = overlaps[0].start
                end = max([x.end for x in overlaps])

            #Negative strand case.
            else:
                start = overlaps[0].start
                end = max([x.end for x in overlaps])

                for overlap in overlaps:
                    if overlap.end == end:
                        name = overlap.value

            #Add all overlapping genes to previously observed overlaps.
            for overlap in overlaps:
                overlapping_genes.add(overlap.value)

        max_dogs.append({'Name':name,'Chrom':dog[1],'Start':start,'End':end,'Strand':dog[4]})

    max_dogs_df = pd.DataFrame(max_dogs)

    return max_dogs_df

'''
Define a function that can find DoG coordinates for pre-screened DoGs given the screening genes, screening coverage, and
DoG window size.
'''
def get_all_dog_coordinates(args):

    screening_genes,screening_coverage_df,window = args

    #Format output for discovering DoGs. Separate genes by strand.
    screening_coverage_copy = screening_coverage_df.copy()
    screening_coverage_copy['ID'] = screening_coverage_copy.apply(lambda x: int(x[0].split('-')[1]),axis=1)
    screening_coverage_copy = pd.merge(screening_coverage_copy,screening_genes[['Name','Strand']],on='Name')
    pos_strand = screening_coverage_copy[screening_coverage_copy.Strand == '+'][['Name','ID']].sort_values(['Name',
                                                                                                            'ID'])
    neg_strand = screening_coverage_copy[screening_coverage_copy.Strand == '-'][['Name','ID']].sort_values(['Name',
                                                                                                            'ID'],
                                                                                                        ascending=False)
    pos_gb = pos_strand.groupby('Name')
    neg_gb = neg_strand.groupby('Name')

    #Get DoG coordinates for all genes.
    dogs = []

    for gene in pos_strand.Name.unique():
        coordinates = get_dog_coordinates(pos_gb.get_group(gene),'+',window)
        dogs.append(coordinates)

    for gene in neg_strand.Name.unique():
        coordinates = get_dog_coordinates(neg_gb.get_group(gene),'-',window)
        dogs.append(coordinates)

    dogs_df = pd.DataFrame(dogs)[['Name','Start','End','Strand']]
    dogs_df = pd.merge(dogs_df,screening_genes[['Name','Chrom']],on='Name')

    #Remove overlaps. Keep DoGs with max coordinates.
    max_dogs_df = max_dogs(dogs_df[['Name','Chrom','Start','End','Strand']])
    max_dogs_df['Score'] = 0
    max_dogs_df = max_dogs_df[['Chrom','Start','End','Name','Score','Strand']]

    return max_dogs_df

'''
Define a function that can output BED files for multiple sets of pre-screened DoGs given a dataframe of screening genes,
a list of screening coverage dataframes, a window, the number of cpus, the tag_directories and an output directory.
'''
def get_multi_dog_beds(screening_genes,screening_coverage_dfs,window,cpu,tag_dirs,out_dir):

    #Get DoGs.
    cmds = []
    for screening_coverage_df in screening_coverage_dfs:
        cmds.append((screening_genes,screening_coverage_df,window))

    pool = Pool(processes=min(cpu,len(screening_coverage_dfs)))
    dog_dfs = pool.map(get_all_dog_coordinates,cmds)
    pool.close()

    #Output BED files.
    for i in range(len(tag_dirs)):
        dog_dfs[i].to_csv(os.path.join(out_dir,tag_dirs[i].split('/')[-1]+'.dogs.bed'),sep='\t',index=False,
                          header=False)

'''
Define a function that can merge all of the DoGs. The rule is that if two or more experiments have the same gene with 
DoG transcript, take the longer transcript.
'''
def merge_dogs(dog_files,out_dir):

    #Load all DoGs into dataframes.
    dog_dfs = []
    for dog_file in dog_files:
        dog_df = pd.read_csv(dog_file,sep='\t',header=None,names=['Chrom','Start','End','Name','Score','Strand'])
        dog_df['Length'] = dog_df.End-dog_df.Start
        dog_dfs.append(dog_df)

    #Merge the dataframes by getting the maximum length for the DoG.
    dog_df = pd.concat(dog_dfs)
    max_length = pd.DataFrame(dog_df.groupby('Name')['Length'].max())
    max_length = max_length.reset_index()
    dog_df = pd.merge(dog_df,max_length,on=['Name','Length']).drop_duplicates()
    del dog_df['Length']

    #Remove overlaps to keep the longest continuous DoG.
    dog_df = max_dogs(dog_df[['Name','Chrom','Start','End','Strand']])
    dog_df['Score'] = 0
    dog_df = dog_df[['Chrom','Start','End','Name','Score','Strand']]

    #Output to BED files.
    dog_df.to_csv(os.path.join(out_dir,'all_dogs.bed'),sep='\t',index=False,header=False)

'''
Define a function that can get the expression for each set of DoGs for each experiment.
'''
def get_dog_exp(tag_dirs,dog_files,stranded,out_dir,cpu):

    cmds = []
    for i in range(len(tag_dirs)):
        cmds.append(([tag_dirs[i]],dog_files[i],stranded,'-raw',out_dir,1))
        cmds.append(([tag_dirs[i]],dog_files[i],stranded,'-fpkm',out_dir,1))

    pool = Pool(processes=min(cpu,len(cmds)))
    pool.map(get_regions_exp,cmds)
    pool.close()

'''
Define a function that can describe the DoG lengths for a single experiment.
'''
def summarize_dog_lens(dogs_bed):

    df = pd.read_csv(dogs_bed,sep='\t',names=['Chrom','Start','End','Name','Score','Strand'])
    df['Length'] = df.End-df.Start
    df = pd.DataFrame(df.describe()['Length'])

    return df

'''
Define a function that can describe the DoG lengths for multiple experiments.
'''
def summarize_dog_lens_all_expts(dogs_beds):

    summary_dfs = []
    for dogs_bed in dogs_beds:
        df = summarize_dog_lens(dogs_bed)
        df.columns = [dogs_bed.split('/')[-1][:-9]]
        summary_dfs.append(df)

    output_df = functools.reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True),summary_dfs)

    return output_df

'''
Define a function that summarizes DoG expression for a single experiment.
'''
def summarize_dog_exp(dog_exp_file):

    df = pd.read_csv(dog_exp_file,sep='\t')
    del df['Length']

    return df.describe()

'''
Define a function that summarizes DoG expression for multiple experiments.
'''
def summarize_dog_exp_all_expts(dog_exp_files):

    summary_dfs = []
    for dog_exp_file in dog_exp_files:
        summary_dfs.append(summarize_dog_exp(dog_exp_file))

    output_df = functools.reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True),summary_dfs)

    return output_df

'''
Create a function that can get the summary for DoGs.
'''
def summarize_all_dogs(all_dogs_bed,dogs_beds,all_dogs_exp,dogs_exp_files,min_dog_len,min_dog_coverage,dog_window):

    summary = f'Summary for DoG finding with minimum length {min_dog_len} bp, minimum coverage of {min_dog_coverage} '+\
              f'FPKM, and screening window of {dog_window} bp\n'

    summary += 'Summary of DoG lengths for all DoGs:\n'+summarize_dog_lens(all_dogs_bed).to_string()

    summary += '\nSummary of DoG lengths for all DoGs for all experiments:\n'+\
               summarize_dog_lens_all_expts(dogs_beds).to_string()

    summary += '\nSummary of expression for all DoGs across all experiments in FPKM\n'+\
               summarize_dog_exp(all_dogs_exp).to_string()

    summary += '\nSummary of expression for DoGs in each experiment in FPKM\n'+\
               summarize_dog_exp_all_expts(dogs_exp_files).to_string()

    return summary