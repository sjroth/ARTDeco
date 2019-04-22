'''
Main script for running ARTDeco. Contains code to run each of the modes.
'''
from .preprocess import parse_gtf,create_stranded_downstream_df,create_stranded_read_in_df,\
    create_unstranded_downstream_df,create_unstranded_read_in_df,make_multi_tag_dirs
from .intergenic import get_multi_gene_exp,get_max_isoform,get_gene_v_intergenic,assign_genes
from .misc import infer_experiments_group,output_inferred_format,get_regions_exp
from .diff_exp_read_in import read_in_diff_exp,assign_read_in_genes
from .get_dogs import get_dog_screening,generate_screening_bed,get_multi_interval_coverage,generate_full_screening_bed,\
    get_multi_dog_beds,merge_dogs,get_dog_exp


import argparse
import os
import sys
import pandas as pd


def main():

    #Make command line interface.
    parser = argparse.ArgumentParser(prog='ARTDeco',
                                     description='Main script for Automatic ReadThrough DEteCtiOn-ARTDeco')
    parser.add_argument('-mode',
                        help='Mode in which to run ARTDeco. Options are preprocess, intergenic, diff_exp_read_in, \
                             get_dogs, and diff_exp_dogs. REQUIRED.',required=True,action='store')
    parser.add_argument('-home-dir',help='Directory in which to run ARTDeco (default is current directory)',
                        action='store',default='.')
    parser.add_argument('-gtf-file',help='GTF file',action='store',default='.')
    parser.add_argument('-cpu',help='Maximum CPUs to use',action='store',type=int,default=1)
    parser.add_argument('-chrom-sizes-file',help='Chromosome sizes file',action='store',default='.')
    parser.add_argument('-read-in-dist',help='Read-in distance. Default is 1 kb.',type=int,default=1000)
    parser.add_argument('-readthrough-dist',help='Readthrough distance. Default is 5 kb.',type=int,default=5000)
    parser.add_argument('-intergenic-min-len',help='Minimum length for intergenic regions. Default is 100 bp.',type=int,
                        default=100)
    parser.add_argument('-intergenic-max-len',help='Maximum length for intergenic regions. Default is 15 kb.',type=int,
                        default=15000)
    parser.add_argument('-read-in-threshold', help='Threshold for considering read-in gene. Default is 0.',type=float,
                        default=0)
    parser.add_argument('-read-in-fpkm',help='Minimum FPKM value for considering a gene. Default is 0.25 FPKM.',
                        type=float,default=0.25)
    parser.add_argument('-overwrite',help='Indicate whether to overwrite existing files',default=False,
                        action='store_true')
    parser.add_argument('-meta-file',help='Meta file',action='store',type=str,default='.')
    parser.add_argument('-comparisons-file',help='Comparisons file',type=str,default='.')
    parser.add_argument('-log2FC', help='Minimum log2 fold change for considering a gene upregulated. Default is 2.',
                        type=float,default=2)
    parser.add_argument('-pval', help='Maximum p-value for considering a gene upregulated. Default is 0.05.',type=float,
                        default=0.05)
    parser.add_argument('-min-dog-len',help='Minimum DoG length. Default is 4 kb.',type=int,default=4000)
    parser.add_argument('-dog-window',help='DoG window size. Default is 500 bp.',type=int,default=500)
    parser.add_argument('-min-dog-coverage',help='Minimum FPKM for DoG discovery. Default is 0.1 FPKM.',type=float,
                        default=0.1)
    args = parser.parse_known_args()[0]

    print(f'Running {args.mode} mode...')

    #Check home directory for BAM files. Store all BAMs and patterns for tag directories.
    no_bams = True
    bam_files = []
    tag_dirs = []
    for f in os.listdir(args.home_dir):
        if f[-4:] == '.bam':
            bam_files.append(f)
            tag_dirs.append(os.path.join(args.home_dir,'preprocess_files',f[:-4]))
            no_bams = False

    #If there are no BAM files, exit.
    if no_bams:
        print('No BAM files... Exiting...')
        sys.exit(1)

    #If the program is running in preprocess, intergenic, or get_dogs mode, infer the format of the BAM files.
    if args.mode in ['preprocess','intergenic','get_dogs']:

        #Check if user-supplied GTF file exists.
        if os.path.isfile(args.gtf_file):
            print('GTF file exists...')
        else:
            print('Invalid GTF file supplied... Exiting...')
            sys.exit(1)

        #Check for existence of a gene annotation directory with all necessary files. Create file if necessary.
        if not os.path.isfile(os.path.join(args.home_dir,'preprocess_files','genes.full.bed')) or \
                not os.path.isfile(os.path.join(args.home_dir,'preprocess_files','genes_condensed.bed')) or \
                not os.path.isfile(os.path.join(args.home_dir,'preprocess_files','gene_to_transcript.txt')) or \
                args.overwrite:

            print('Incomplete gene annotation files or overwrite specified... Will create gene annotations...')

            #Create directory if it does not already exist.
            if not os.path.isdir(os.path.join(args.home_dir,'preprocess_files')):
                os.mkdir(os.path.join(args.home_dir,'preprocess_files'))

            #Create the annotation files.
            parse_gtf(gtf_file=args.gtf_file,home_dir=args.home_dir)

        else:
            print('Gene annotation files exist...')

        print('Inferring BAM file formats...')
        formats = infer_experiments_group([os.path.join(args.home_dir,bam_file) for bam_file in bam_files],
                                          os.path.join(args.home_dir,'preprocess_files','genes.full.bed'),
                                          min(args.cpu,len(bam_files)))

        #Check that all of the BAM files are of the same format.
        if len(set(x[1] for x in formats)) == 1 and len(set(x[2] for x in formats)) == 1 \
                and len(set(x[3] for x in formats)) == 1:

            pe = formats[0][1]
            stranded = formats[0][2]
            flip = formats[0][3]

            out_str = 'All BAM files are ' + output_inferred_format(formats[0])

            print(out_str)

        else:
            print('Error... One or more files do not match in inferred format... Exiting...')

            for f in formats:
                out_str = f'BAM file {f[0]} inferred as '+output_inferred_format(f)
                print(out_str)
                sys.exit(1)

    #If the program is running in intergenic or get_dogs mode, check for existence of tag directories.
    if args.mode in ['intergenic','get_dogs']:

        for tag_dir in tag_dirs:

            if not os.path.isdir(tag_dir):
                print('Tag directories missing... Please run preprocessing mode... Exiting...')
                sys.exit(1)

        print('All tag directories exist...')

    #If the program is running in diff_exp_read_in or diff_exp_dogs mode, check for properly formatted meta and
    #comparisons files.
    if args.mode in ['diff_exp_read_in','diff_exp_dogs']:

        #Import necessary packages.
        from .DESeq2 import reformat_meta,reformat_comparisons,generate_comparisons,load_deseq_dataset,run_deseq,\
            deseq_results

        #Check if reformatted files already exist. If they do not, generate them.
        if os.path.isfile(os.path.join(args.home_dir,'preprocess_files','meta.reformatted.txt')) and \
                os.path.isfile(os.path.join(args.home_dir,'preprocess_files','comparisons.reformatted.txt')) and \
                not args.overwrite:
            print('Reformatted meta and comparisons files exist...')

        else:

            #Check meta file.
            if os.path.isfile(args.meta_file):

                meta = open(args.meta_file).readlines()
                meta_format = True
                groups = []
                i = 0
                while i < len(meta) and meta_format:

                    line = meta[i].strip().split('\t')

                    #If the length of the split line is different than 2, the meta file isn't properly formatted.
                    if len(line) != 2:
                        meta_format = False

                    #Otherwise, ensure that the first line is the proper format.
                    else:
                        if i == 0:
                            if line[0] != 'Experiment' or line[1] != 'Group':
                                meta_format = False
                        else:
                            if line[1] not in groups:
                                groups.append(line[1])

                    i += 1

                if meta_format:
                    print('Meta file properly formatted... Generating reformatted meta...')
                    reformat_meta(args.meta_file,os.path.join(args.home_dir,'preprocess_files'))
                else:
                    print('Meta file not properly formatted... Exiting...')
                    sys.exit(1)

            else:
                print('Meta file does not exist... Exiting...')
                sys.exit(1)

            #Check comparisons file.
            if os.path.isfile(args.comparisons_file):

                print('Comparisons file exists...')

                #Check format.
                comparisons = [line.strip().split('\t') for line in open(args.comparisons_file).readlines()]
                comparisons_lens = [len(line) for line in comparisons]

                #Check if lines are tab-separated formatted.
                if len(set(comparisons_lens)) == 1 and len(comparisons[0]) == 2:

                    membership = [(line[0] in groups and line[1] in groups) for line in comparisons]

                    if len(set(membership)) == 1 and membership[0]:
                        comparisons_format = True
                    else:
                        comparisons_format = False
                else:
                    comparisons_format = False

                #If the file is properly formatted, reformat it. Otherwise, generate an all-by-all file.
                if comparisons_format:
                    print('Comparisons file properly formatted... Generating reformatted comparisons...')
                    reformat_comparisons(args.comparisons_file,os.path.join(args.home_dir,'preprocess_files'))
                else:
                    print('Comparisons file not properly formatted... Generating all-by-all comparisons file...')
                    generate_comparisons(os.path.join(args.home_dir,'preprocess_files','meta.reformatted.txt'),
                                         os.path.join(args.home_dir,'preprocess_files'))

            else:
                print('Comparison file does not exist or not provided... Generating comparisons file...')
                generate_comparisons(os.path.join(args.home_dir,'preprocess_files','meta.reformatted.txt'),
                                     os.path.join(args.home_dir,'preprocess_files'))

        comparisons = [line.strip().split('\t') for line in
                       open(os.path.join(args.home_dir,'preprocess_files','comparisons.reformatted.txt')).readlines()]

    #If the program is running in preprocess or get_dogs mode, check if the chromosome sizes file exists.
    #Load chromosome sizes file.
    if args.mode in ['preprocess','get_dogs']:
        if os.path.isfile(args.chrom_sizes_file):
            chrom_sizes_file = open(args.chrom_sizes_file)
            chrom_sizes = {}
            for line in chrom_sizes_file.readlines():
                line = line.strip().split('\t')
                chrom_sizes[line[0]] = int(line[1])
        else:
            print('Chromosome sizes file does not exist... Exiting...')
            sys.exit(1)

    #Run preprocess mode.
    if args.mode == 'preprocess':

        #Create read-in and downstream BED files.
        if not os.path.isfile(os.path.join(args.home_dir,'preprocess_files','read_in.bed')) or \
                not os.path.isfile(os.path.join(args.home_dir,'preprocess_files','readthrough.bed')) or args.overwrite:

            print('Creating intergenic BED files...')

            #Load genes.
            genes = pd.read_csv(os.path.join(args.home_dir,'preprocess_files/genes_condensed.bed'),sep='\t',header=None,
                                names=['Chrom','Start','Stop','Name','Score','Strand'])
            del genes['Score']
            genes = genes[['Name','Chrom','Strand','Start','Stop']]

            #Create intergenic dataframes.
            if stranded:
                read_in_df = create_stranded_read_in_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                        min_len=args.intergenic_min_len,upstream_dist=args.read_in_dist)
                readthrough_df = create_stranded_downstream_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                               min_len=args.intergenic_min_len,
                                                               downstream_dist=args.readthrough_dist)
            else:
                read_in_df = create_unstranded_read_in_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                          min_len=args.intergenic_min_len,
                                                          upstream_dist=args.read_in_dist)
                readthrough_df = create_unstranded_downstream_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                                 min_len=args.intergenic_min_len,
                                                                 downstream_dist=args.readthrough_dist)

            #Format and create BED files.
            read_in_df['Score'] = 0
            read_in_df = read_in_df[['Chrom','Start','Stop','Name','Score','Strand']]
            read_in_df['Start'] = read_in_df['Start'].astype(int)
            read_in_df['Stop'] = read_in_df['Stop'].astype(int)

            readthrough_df['Score'] = 0
            readthrough_df = readthrough_df[['Chrom','Start','Stop','Name','Score','Strand']]
            readthrough_df['Start'] = readthrough_df['Start'].astype(int)
            readthrough_df['Stop'] = readthrough_df['Stop'].astype(int)

            read_in_df.to_csv(os.path.join(args.home_dir,'preprocess_files','read_in.bed'),sep='\t',header=False,
                              index=False)
            readthrough_df.to_csv(os.path.join(args.home_dir,'preprocess_files','readthrough.bed'),sep='\t',
                                  header=False,index=False)

        else:
            print('Intergenic files already exist...')

        #Create tag directories if they do not exist.
        tag_dir_cmds = []
        for i in range(len(tag_dirs)):
            if not os.path.isdir(tag_dirs[i]) or args.overwrite:
                tag_dir_cmds.append(os.path.join(args.home_dir,bam_files[i]))

        if tag_dir_cmds:
            print('Creating tag directories...')
            make_multi_tag_dirs(tag_dir_cmds,os.path.join(args.home_dir,'preprocess_files'),flip,pe,stranded,
                                min(len(bam_files),args.cpu))
        else:
            print('Tag directories already exist...')

    #Run intergenic mode.
    if args.mode == 'intergenic':

        #Check if quantification and intergenic directories exist. Create them if they do not.
        if not os.path.isdir(os.path.join(args.home_dir,'quantification')):
            os.mkdir(os.path.join(os.path.join(args.home_dir,'quantification')))
        if not os.path.isdir(os.path.join(args.home_dir,'intergenic')):
            os.mkdir(os.path.join(os.path.join(args.home_dir,'intergenic')))

        #Check if gene expression and maximum isoform files exist. If one of them does not, create all three.
        if os.path.isfile(os.path.join(args.home_dir,'quantification','gene.exp.raw.txt')) and \
                os.path.isfile(os.path.join(args.home_dir,'quantification','gene.exp.fpkm.txt')) and \
                os.path.isfile(os.path.join(args.home_dir,'quantification','max_isoform.txt')) and not args.overwrite:
            print('Gene expression files exist...')
        else:
            print('Generating gene expression and maximum isoform files...')
            get_multi_gene_exp(tag_dirs,args.gtf_file,stranded,os.path.join(args.home_dir,'quantification'),args.cpu)
            get_max_isoform(os.path.join(args.home_dir,'quantification','gene.exp.fpkm.txt'),
                            os.path.join(args.home_dir,'preprocess_files','gene_to_transcript.txt'),
                            os.path.join(args.home_dir,'quantification'))

        #Check if read-in and readthrough expression files exist. If one of them does not, create both of them.
        if os.path.isfile(os.path.join(args.home_dir,'quantification','read_in.exp.txt')) and \
                os.path.isfile(os.path.join(args.home_dir,'quantification','readthrough.exp.txt')) and not args.overwrite:
            print('Intergenic expression files exist...')
        else:
            print('Generating intergenic expression files...')
            get_regions_exp((tag_dirs,os.path.join(args.home_dir,'preprocess_files','read_in.bed'),stranded,'-raw',
                             os.path.join(args.home_dir,'quantification'),min(len(tag_dirs),args.cpu)))
            get_regions_exp((tag_dirs,os.path.join(args.home_dir,'preprocess_files', 'readthrough.bed'),stranded,'-raw',
                             os.path.join(args.home_dir,'quantification'),min(len(tag_dirs),args.cpu)))

        #Check if intergenic vs. gene expression files exist. IF one of them does not, create both of them.
        if os.path.isfile(os.path.join(args.home_dir,'intergenic','read_in.txt')) and \
                os.path.isfile(os.path.join(args.home_dir,'intergenic','readthrough.txt')) and \
                os.path.isfile(os.path.join(args.home_dir,'intergenic','read_in_assignments.txt')) and \
                not args.overwrite:
            print('Intergenic vs. gene expression files and read-in gene assignments exist...')
        else:
            print(f'Generating intergenic vs. expression files and read-in gene assignments... Read-in level '+\
                  f'threshold is {args.read_in_threshold} and read-in FPKM threshold is {args.read_in_fpkm}...')

            get_gene_v_intergenic(os.path.join(args.home_dir,'quantification','gene.exp.raw.txt'),
                                  os.path.join(args.home_dir,'quantification','gene.exp.fpkm.txt'),
                                  os.path.join(args.home_dir,'quantification','max_isoform.txt'),
                                  os.path.join(args.home_dir,'quantification','read_in.raw.txt'),'Read-In',
                                  os.path.join(args.home_dir,'intergenic','read_in.txt'))

            get_gene_v_intergenic(os.path.join(args.home_dir,'quantification','gene.exp.raw.txt'),
                                  os.path.join(args.home_dir,'quantification','gene.exp.fpkm.txt'),
                                  os.path.join(args.home_dir,'quantification','max_isoform.txt'),
                                  os.path.join(args.home_dir,'quantification','readthrough.raw.txt'),'Readthrough',
                                  os.path.join(args.home_dir,'intergenic','readthrough.txt'))

            assign_genes(os.path.join(args.home_dir,'intergenic','read_in.txt'),args.read_in_threshold,
                         args.read_in_fpkm,os.path.join(args.home_dir,'intergenic','read_in_assignments.txt'))

    #Run diff_exp_read_in mode.
    if args.mode == 'diff_exp_read_in':

        #Run DESeq2.
        print('Check for DESeq2 output...')
        all_comparisons = True
        for condition1,condition2 in comparisons:
            if not os.path.isfile(os.path.join(args.home_dir,'diff_exp',f'{condition1}-{condition2}-results.txt')):
                all_comparisons = False

        if all_comparisons and not args.overwrite:
            print('DESeq2 results exist...')
        else:

            print('Running DESeq2')
            dds = load_deseq_dataset(os.path.join(args.home_dir,'quantification','gene.exp.raw.txt'),
                                     os.path.join(args.home_dir,'preprocess_files','meta.reformatted.txt'))
            dds_results = run_deseq(dds)

            #Create differential expression directory.
            if not os.path.isdir(os.path.join(args.home_dir,'diff_exp')):
                os.mkdir(os.path.join(args.home_dir,'diff_exp'))

            #Output results.
            print('Output DESeq2 results...')
            for condition1,condition2 in comparisons:
                deseq_results(dds_results,condition1,condition2,os.path.join(args.home_dir,'diff_exp'))

        #Create differential expression with read-in information directory.
        print('Check for differential expression joined with read-in information as well as read-in gene '+\
              'assignments...')
        all_read_in = True
        for condition1,condition2 in comparisons:
            if not os.path.isfile(os.path.join(args.home_dir,'diff_exp_read_in',
                                               f'{condition1}-{condition2}-read_in.txt')) or \
                    not os.path.isfile(os.path.join(args.home_dir,'diff_exp_read_in',
                                                    f'{condition1}-{condition2}-read_in_assignment.txt')):
                all_read_in = False

        if all_read_in and not args.overwrite:
            print('Differential expression with read-in information and read-in gene assignments exist...')
        else:

            if not os.path.isdir(os.path.join(args.home_dir,'diff_exp_read_in')):
                os.mkdir(os.path.join(args.home_dir,'diff_exp_read_in'))

            #Join differential expression information to read-in information. Use this to infer read-in genes.
            print('Combining differential expression results and read-in information... Inferring read-in genes for '+\
                  f'upregulated genes with log2 fold change > {args.log2FC}, p-value < {args.pval}, and FPKM > '+\
                  f'{args.read_in_fpkm}... Read-in level threshold is {args.read_in_threshold}...')

            for condition1,condition2 in comparisons:
                read_in_diff_exp(os.path.join(args.home_dir,'intergenic/read_in.txt'),
                                 os.path.join(args.home_dir,'preprocess_files/meta.reformatted.txt'),
                                 os.path.join(args.home_dir,'diff_exp', f'{condition1}-{condition2}-results.txt'),
                                 os.path.join(args.home_dir,'diff_exp_read_in'))

                assign_read_in_genes(os.path.join(args.home_dir,'diff_exp_read_in',
                                                  f'{condition1}-{condition2}-read_in.txt'),
                                     args.log2FC,args.pval,args.read_in_fpkm,args.read_in_threshold,
                                     os.path.join(args.home_dir,'diff_exp_read_in'))

    #Run get_dogs mode.
    if args.mode == 'get_dogs':

        #Create DoG directory if one does not exist.
        if not os.path.isdir(os.path.join(args.home_dir,'dogs')):
           os.mkdir(os.path.join(args.home_dir,'dogs'))

        dog_files = []
        all_bed_exist = True
        all_exp_exist = True
        for tag_dir in tag_dirs:

            dog = os.path.join(args.home_dir,'dogs',tag_dir.split('/')[-1])+'.dogs.'

            dog_files.append(dog+'bed')

            if not os.path.isfile(dog+'bed'):
                all_bed_exist = False

            if not os.path.isfile(dog+'raw.txt') or not os.path.isfile(dog+'fpkm.txt'):
                all_exp_exist = False

        if not os.path.isfile(os.path.join(args.home_dir,'dogs','all_dogs.bed')):
            all_bed_exist = False

        if not os.path.isfile(os.path.join(args.home_dir,'dogs','all_dogs.raw.txt')) or \
                not os.path.isfile(os.path.join(args.home_dir,'dogs','all_dogs.fpkm.txt')):
            all_exp_exist = False

        if not all_bed_exist or args.overwrite:

            print('Finding DoGs...')

            print(f'Get genes with potential DoGs with minimum length of {args.min_dog_len} bp...')
            screening_genes = get_dog_screening(os.path.join(args.home_dir,'preprocess_files','genes_condensed.bed'),
                                                args.min_dog_len)

            print(f'Generate initial screening BED file for DoGs with minimum length {args.min_dog_len} bp and window'+\
                  f' size {args.dog_window} bp...')
            generate_screening_bed(screening_genes,args.min_dog_len,args.dog_window,args.home_dir)

            print(f'Initial screening coverage for DoGs with minimum length of {args.min_dog_len} bp...')

            #Screen for coverage threshold.
            screening_coverage_dfs = get_multi_interval_coverage(tag_dirs,args.home_dir,
                                                                 [os.path.join(args.home_dir,'intervals.bed')],
                                                                 args.min_dog_coverage,stranded,min(args.cpu,len(tag_dirs)))

            #Find genes that pass threshold for minimum length.
            screening_genes_dfs = {}
            for i in range(len(screening_coverage_dfs)):
                new_df = screening_coverage_dfs[i].copy()
                new_df = new_df.groupby('Name').count()
                new_df = new_df[new_df.ID == args.min_dog_len/args.dog_window]
                screening_genes_dfs[tag_dirs[i]] = screening_genes[screening_genes.Name.isin(new_df.index)].copy()

            #Remove screening BED files.
            os.remove(os.path.join(args.home_dir,'intervals.bed'))

            print('Generate screening BED file for pre-screened DoGs...')

            generate_full_screening_bed(tag_dirs, os.path.join(args.home_dir,'preprocess_files','genes_condensed.bed'),
                                        screening_genes_dfs,
                                        os.path.join(args.home_dir,'intergenic','read_in_assignments.txt'),
                                        args.chrom_sizes_file,args.dog_window,args.cpu,args.home_dir)

            print('Screening coverage for pre-screened DoGs...')

            #Screen for coverage threshold.
            screening_coverage_dfs = get_multi_interval_coverage(tag_dirs,args.home_dir,
                                                                 [os.path.join(args.home_dir,
                                                                               tag_dir.split('/')[-1]+'.bed') for
                                                                  tag_dir in tag_dirs],args.min_dog_coverage,stranded,
                                                                 min(args.cpu,len(tag_dirs)))

            #Remove screening BED files.
            for tag_dir in tag_dirs:
                expt_name = tag_dir.split('/')[-1]
                os.remove(os.path.join(args.home_dir,f'{expt_name}.bed'))

            print('Discovering DoG coordinates for pre-screened DoGs and output BED files...')
            get_multi_dog_beds(screening_genes,screening_coverage_dfs,args.dog_window,args.cpu,tag_dirs,
                               os.path.join(args.home_dir,'dogs'))

            print('Merge DoGs into a single annotation...')
            merge_dogs(dog_files,os.path.join(args.home_dir,'dogs'))

        else:
            print('All DoG annotations exist...')

        if not all_exp_exist or args.overwrite:

            print('Quantifying DoGs...')

            print('Quantify DoGs for each experiment...')
            get_dog_exp(tag_dirs,dog_files,stranded,os.path.join(args.home_dir,'dogs'),args.cpu)

            print('Quantify all DoGs for all experiments...')
            get_regions_exp((tag_dirs,os.path.join(args.home_dir,'dogs','all_dogs.bed'),stranded,'-raw',
                             os.path.join(args.home_dir,'dogs'),min(args.cpu,len(tag_dirs))))
            get_regions_exp((tag_dirs, os.path.join(args.home_dir,'dogs','all_dogs.bed'),stranded,'-fpkm',
                             os.path.join(args.home_dir,'dogs'),min(args.cpu,len(tag_dirs))))
        else:
            print('DoG quantification files exist...')

    #Run diff_exp_dogs mode.
    if args.mode == 'diff_exp_dogs':

        #Run DESeq2.
        print('Check for DESeq2 output...')
        all_comparisons = True
        for condition1, condition2 in comparisons:
            if not os.path.isfile(os.path.join(args.home_dir,'diff_exp_dogs',f'{condition1}-{condition2}-results.txt')):
                all_comparisons = False

        if all_comparisons and not args.overwrite:
            print('DESeq2 results exist...')
        else:
            print('Running DESeq2')
            dds = load_deseq_dataset(os.path.join(args.home_dir,'dogs','all_dogs.raw.txt'),
                                     os.path.join(args.home_dir,'preprocess_files','meta.reformatted.txt'))
            dds_results = run_deseq(dds)

            #Create differential expression directory.
            if not os.path.isdir(os.path.join(args.home_dir,'diff_exp_dogs')):
                os.mkdir(os.path.join(args.home_dir,'diff_exp_dogs'))

            #Output results.
            print('Output DESeq2 results...')
            for condition1,condition2 in comparisons:
                deseq_results(dds_results,condition1,condition2,os.path.join(args.home_dir,'diff_exp_dogs'))

if __name__ == '__main__':
    main()