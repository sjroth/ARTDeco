'''
Module that contains assorted utils that do not fall under a single usage category. Not including DESeq2 stuff as it
takes time to load.
'''
import networkx as nx
import subprocess
from multiprocessing import Pool
import pandas as pd
import os

'''
Define an object that manages file dependencies for the ARTDeco output. 
'''
class ARTDecoDir:

    #Initialization method. Creates a dependency graph.
    def __init__(self,bam_files_dir,home_dir):

        self.bam_files_dir = bam_files_dir
        self.home_dir = home_dir

        #Create graph.
        self.dependency = nx.DiGraph()

        #ARTDeco directories.
        self.preprocess_dir = os.path.join(home_dir,'preprocess_files')
        self.quantification_dir = os.path.join(home_dir,'quantification')
        self.readthrough_dir = os.path.join(home_dir,'readthrough')
        self.dogs_dir = os.path.join(home_dir,'dogs')
        self.diff_exp_dir = os.path.join(home_dir,'diff_exp')
        self.diff_exp_read_in_dir = os.path.join(home_dir,'diff_exp_read_in')
        self.diff_exp_dogs_dir = os.path.join(home_dir,'diff_exp_dogs')
        self.summary_dir = os.path.join(home_dir,'summary_files')

        #Add nodes.
        self.dependency.add_nodes_from([self.preprocess_dir,self.quantification_dir,self.readthrough_dir,self.dogs_dir,
                                        self.diff_exp_dir,self.diff_exp_read_in_dir,self.diff_exp_dogs_dir])

        #Grab BAM files.
        self.bam_files = []
        self.tag_dirs = set()
        self.all_tag_dirs = []
        self.tag_dir_to_bam = {}
        for f in os.listdir(bam_files_dir):
            if f[-4:] == '.bam':
                bam_file = os.path.join(bam_files_dir,f)
                tag_dir = os.path.join(self.preprocess_dir,f[:-4])

                self.bam_files.append(bam_file)
                self.tag_dirs.add(tag_dir)
                self.all_tag_dirs.append(tag_dir)
                self.tag_dir_to_bam[tag_dir] = bam_file

        self.all_tag_dirs.sort()

        #Preprocessing files.
        self.comparisons_file = os.path.join(self.preprocess_dir,'comparisons.reformatted.txt')
        self.genes_condensed = os.path.join(self.preprocess_dir,'genes_condensed.bed')
        self.genes_full = os.path.join(self.preprocess_dir,'genes.full.bed')
        self.gene_to_transcript = os.path.join(self.preprocess_dir,'gene_to_transcript.txt')
        self.gene_types = os.path.join(self.preprocess_dir,'gene_types.txt')
        self.meta_file = os.path.join(self.preprocess_dir,'meta.reformatted.txt')
        self.read_in_bed = os.path.join(self.preprocess_dir,'read_in.bed')
        self.readthrough_bed = os.path.join(self.preprocess_dir,'readthrough.bed')

        self.preprocessing_files = {self.genes_condensed,self.genes_full,self.gene_to_transcript,self.gene_types,
                                    self.read_in_bed,self.readthrough_bed} | self.tag_dirs

        #Add nodes.
        self.dependency.add_nodes_from(self.preprocessing_files)

        #Add edges.
        self.dependency.add_edges_from([(self.preprocess_dir,self.comparisons_file),
                                        (self.preprocess_dir,self.genes_full),(self.preprocess_dir,self.meta_file),
                                        (self.genes_full,self.genes_condensed),
                                        (self.genes_full,self.gene_types),
                                        (self.genes_condensed,self.gene_to_transcript),
                                        (self.genes_condensed,self.read_in_bed),
                                        (self.genes_condensed,self.readthrough_bed)]+
                                       [(self.preprocess_dir,tag_dir) for tag_dir in self.tag_dirs])

        #Quantification files.
        self.gene_fpkm = os.path.join(self.quantification_dir,'gene.exp.fpkm.txt')
        self.gene_raw = os.path.join(self.quantification_dir,'gene.exp.raw.txt')
        self.max_isoform = os.path.join(self.quantification_dir,'max_isoform.txt')
        self.read_in_exp = os.path.join(self.quantification_dir,'read_in.raw.txt')
        self.readthrough_exp = os.path.join(self.quantification_dir,'readthrough.raw.txt')

        self.quantification_files = {self.gene_fpkm,self.gene_raw,self.max_isoform,self.read_in_exp,
                                     self.readthrough_exp}

        #Add nodes.
        self.dependency.add_nodes_from(self.quantification_files)

        #Add edges.
        self.dependency.add_edges_from([(self.quantification_dir,self.gene_fpkm),
                                        (self.quantification_dir,self.gene_raw),
                                        (self.quantification_dir,self.read_in_exp),
                                        (self.quantification_dir,self.readthrough_exp)]+
                                       [(tag_dir,self.gene_fpkm) for tag_dir in self.tag_dirs]+
                                       [(tag_dir,self.gene_raw) for tag_dir in self.tag_dirs]+
                                       [(tag_dir,self.read_in_exp) for tag_dir in self.tag_dirs]+
                                       [(tag_dir,self.readthrough_exp) for tag_dir in self.tag_dirs]+
                                       [(self.gene_fpkm,self.max_isoform),(self.gene_to_transcript,self.max_isoform),
                                        (self.read_in_bed,self.read_in_exp),(self.readthrough_bed,self.readthrough_exp),
                                        (self.read_in_bed,self.read_in_exp)])

        #Readthrough files.
        self.corrected_exp = os.path.join(self.readthrough_dir, 'corrected_exp.txt')
        self.read_in_assignments = os.path.join(self.readthrough_dir,'read_in_assignments.txt')
        self.read_in_levels = os.path.join(self.readthrough_dir,'read_in.txt')
        self.readthrough_levels = os.path.join(self.readthrough_dir,'readthrough.txt')

        self.readthrough_files = {self.read_in_assignments,self.read_in_levels,self.readthrough_levels,
                                  self.corrected_exp}

        #Add nodes.
        self.dependency.add_nodes_from(self.readthrough_files)

        #Add edges.
        self.dependency.add_edges_from([(self.readthrough_dir,self.corrected_exp),
                                        (self.readthrough_dir,self.read_in_levels),
                                        (self.readthrough_dir,self.readthrough_levels),
                                        (self.gene_fpkm,self.read_in_levels),(self.gene_raw,self.read_in_levels),
                                        (self.max_isoform,self.read_in_levels),(self.read_in_exp,self.read_in_levels),
                                        (self.gene_fpkm,self.readthrough_levels),
                                        (self.gene_raw,self.readthrough_levels),
                                        (self.max_isoform,self.readthrough_levels),
                                        (self.readthrough_exp,self.readthrough_levels),
                                        (self.read_in_levels,self.corrected_exp),
                                        (self.read_in_levels,self.read_in_assignments),
                                        (self.gene_types,self.read_in_assignments)])

        #DoGs files.
        self.all_dogs_bed = os.path.join(self.dogs_dir,'all_dogs.bed')
        self.all_dogs_fpkm = os.path.join(self.dogs_dir,'all_dogs.fpkm.txt')
        self.all_dogs_raw = os.path.join(self.dogs_dir,'all_dogs.raw.txt')

        self.dogs_beds = set()
        self.all_dogs_beds = []
        self.dogs_bed_to_tagdir = {}
        self.dogs_raw = set()
        self.dogs_fpkm = set()
        self.all_dogs_fpkm_expts = []
        for tag_dir in self.tag_dirs:
            dog = os.path.join(home_dir,'dogs',tag_dir.split('/')[-1])+'.dogs.'

            self.dogs_beds.add(dog+'bed')
            self.all_dogs_beds.append(dog+'bed')
            self.dogs_bed_to_tagdir[dog+'bed'] = tag_dir
            self.dogs_raw.add(dog+'raw.txt')
            self.dogs_fpkm.add(dog+'fpkm.txt')
            self.all_dogs_fpkm_expts.append(dog+'fpkm.txt')

            self.dependency.add_nodes_from([dog+'bed',dog+'raw.txt',dog+'fpkm.txt'])
            self.dependency.add_edges_from([(self.dogs_dir,dog+'bed'),(self.genes_condensed,dog+'bed'),
                                            (tag_dir,dog+'bed'),(self.read_in_assignments,dog+'bed'),
                                            (tag_dir,dog+'raw.txt'),(tag_dir,dog+'fpkm.txt'),(dog+'bed',dog+'raw.txt'),
                                            (dog+'bed',dog+'fpkm.txt')])

        self.all_dogs_beds.sort()

        self.dogs_files = {self.all_dogs_bed,self.all_dogs_fpkm,self.all_dogs_raw} | self.dogs_beds | self.dogs_raw |\
                          self.dogs_fpkm

        #Add nodes.
        self.dependency.add_nodes_from([self.all_dogs_bed,self.all_dogs_fpkm,self.all_dogs_raw])

        #Add edges.
        self.dependency.add_edges_from([(dog,self.all_dogs_bed) for dog in self.dogs_beds]+
                                       [(tag_dir,self.all_dogs_fpkm) for tag_dir in self.tag_dirs]+
                                       [(tag_dir,self.all_dogs_raw) for tag_dir in self.tag_dirs]+
                                       [(self.all_dogs_bed,self.all_dogs_fpkm),(self.all_dogs_bed,self.all_dogs_raw)])

        #Assorted lists for prerequisite files.
        self.gtf_needed = {self.genes_condensed,self.genes_full,self.gene_to_transcript,self.gene_types,self.gene_fpkm,
                           self.gene_raw}
        self.format_needed = {self.read_in_bed,self.readthrough_bed,self.gene_fpkm,self.gene_raw,self.read_in_exp,
                              self.readthrough_exp,self.all_dogs_fpkm,self.all_dogs_raw} | self.tag_dirs | \
                             self.dogs_beds | self.dogs_fpkm | self.dogs_fpkm
        self.chrom_sizes_needed = {self.read_in_bed,self.readthrough_bed} | self.dogs_beds

    #Define a function that can set up differential expression outputs if prompted. Assumes comparisons file has
    #been generated.
    def set_diff_exp_output(self):

        self.__comparisons = [line.strip().split('\t') for line in open(self.comparisons_file).readlines()]

        #Differential expression output files.
        self.diff_exp_files = set()
        self.diff_exp_read_in_files = set()
        self.diff_exp_dogs_files = set()
        for condition1,condition2 in self.__comparisons:
            diff_exp = os.path.join(self.diff_exp_dir,f'{condition1}-{condition2}-results.txt')
            diff_exp_read_in = os.path.join(self.diff_exp_read_in_dir,f'{condition1}-{condition2}-read_in.txt')
            assignments = os.path.join(self.diff_exp_read_in_dir,f'{condition1}-{condition2}-read_in_assignment.txt')
            diff_exp_dogs = os.path.join(self.diff_exp_dogs_dir,f'{condition1}-{condition2}-results.txt')

            self.diff_exp_files.add(diff_exp)
            self.diff_exp_read_in_files |= {diff_exp_read_in,assignments}
            self.diff_exp_dogs_files.add(diff_exp_dogs)

            #Add nodes.
            self.dependency.add_nodes_from([diff_exp,diff_exp_read_in,assignments,diff_exp_dogs])

            #Add edges.
            self.dependency.add_edges_from([(self.diff_exp_dir,diff_exp),(self.gene_raw,diff_exp),
                                            (self.diff_exp_read_in_dir,diff_exp_read_in),
                                            (self.read_in_levels,diff_exp_read_in),(diff_exp,diff_exp_read_in),
                                            (diff_exp_read_in,assignments),(self.diff_exp_dogs_dir,diff_exp_dogs),
                                            (self.all_dogs_raw,diff_exp_dogs)])

    #Get directories/files to generate based upon mode.
    def get_files(self,mode,meta_provided,overwrite):

        flatten = lambda l: [item for sublist in l for item in sublist]

        #If mode provided, specify output files.
        if mode:
            if mode == 'preprocess':
                end_files = self.preprocessing_files
                if not meta_provided:
                    end_files.remove(self.comparisons_file)
                    end_files.remove(self.meta_file)
            elif mode == 'readthrough':
                end_files = self.readthrough_files
            elif mode == 'get_dogs':
                end_files = self.dogs_files
            elif mode == 'diff_exp_read_in':
                end_files = self.diff_exp_read_in_files
            else:
                end_files = self.diff_exp_dogs_files

        #No mode but meta provided.
        elif meta_provided:
            end_files = self.preprocessing_files | self.readthrough_files | self.dogs_files | \
                        self.diff_exp_read_in_files | self.diff_exp_dogs_files

        #No mode with no meta.
        else:
            end_files = self.preprocessing_files | self.readthrough_files | self.dogs_files

        #Go through paths and get non-existent out files.
        out_files = set()
        for f in end_files:
            paths = nx.single_target_shortest_path(self.dependency,f)
            for out_file in set(flatten(paths.values())):
                if overwrite:
                    if not os.path.isdir(out_file) or out_file in self.tag_dirs:
                        out_files.add(out_file)
                else:
                    if not os.path.isfile(out_file) and not os.path.isdir(out_file):
                        out_files.add(out_file)

        return out_files

    #Update file directory lists to reflect what needs to be generated.
    def update_dir_lists(self,out_files):
        self.tag_dirs &= out_files
        self.preprocessing_files &= out_files
        self.quantification_files &= out_files
        self.readthrough_files &= out_files
        self.dogs_files &= out_files
        self.dogs_beds &= out_files
        try:
            self.diff_exp_files &= out_files
            self.diff_exp_read_in_files &= out_files
            self.diff_exp_dogs_files &= out_files
        except:
            pass
        self.gtf_needed &= out_files
        self.format_needed &= out_files
        self.chrom_sizes_needed &= out_files

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
Define a function that can count the total number of reads given a flag for a BAM file.
'''
def count_reads(args):

    bam_file,flag = args

    if flag:
        flag_cmd = [flag]
    else:
        flag_cmd = []

    p = subprocess.Popen(['samtools','view','-c']+flag_cmd+[bam_file],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output,err = p.communicate()

    return int(output.strip())

'''
Define a function that can summarize file formats.
'''
def summarize_bam_files(bam_files,cpu,pe,stranded,flip):

    #Get read counts.
    cmds = []
    for bam_file in bam_files:
        cmds += [(bam_file,''),(bam_file,'-F 4')]
    pool = Pool(processes=min(cpu,len(cmds)))
    read_counts = pool.map(count_reads,cmds)
    pool.close()

    #Create file summary text.
    out_text = f'{len(bam_files)} Experiments\n'

    if pe:
        pe_str = 'Paired-End'
    else:
        pe_str = 'Single-End'
    if stranded:
        stranded_str = 'Strand-Specific'
        if flip:
            orientation_str = 'Reverse-strand oriented'
        else:
            orientation_str = 'Forward-strand oriented'
    else:
        stranded_str = 'Unstranded'
        orientation_str = ''

    out_text += f'Files are {pe_str}, {stranded_str}'
    if orientation_str:
        out_text += f', {orientation_str}\n'
    else:
        out_text += '\n'

    out_lst = []
    for i in range(len(bam_files)):
        out_lst.append(
            {'Experiment':bam_files[i],'Total Reads':read_counts[2*i],'Mapped Reads':read_counts[2*i+1]})

    out_df = pd.DataFrame(out_lst)
    out_df = out_df[['Experiment','Total Reads','Mapped Reads']]

    return out_text+'\n'.join(out_df.to_string(index=False).split('\n'))

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
            new_name = new_name.replace(' ','_')
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