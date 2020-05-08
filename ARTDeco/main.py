'''
Main script for running ARTDeco. Contains code to run each of the modes.
'''
from .misc import ARTDecoDir,infer_experiments_group,output_inferred_format,summarize_bam_files,get_regions_exp
from .preprocess import parse_gtf,create_stranded_downstream_df,create_stranded_read_in_df,\
    create_unstranded_downstream_df,create_unstranded_read_in_df,make_multi_tag_dirs
from .readthrough import get_multi_gene_exp,get_max_isoform,get_gene_v_intergenic,deconvolute_exp,assign_genes,\
    summarize_readthrough_stats,summarize_read_in_assignments
from .diff_exp_read_in import read_in_diff_exp,assign_read_in_genes,summarize_diff_exp_read_in_assignments
from .get_dogs import get_dog_screening,generate_screening_bed,get_multi_interval_coverage,generate_full_screening_bed,\
    get_multi_dog_beds,merge_dogs,get_dog_exp,summarize_all_dogs


import argparse
import os
import sys
import pandas as pd


def main():

    #Make command line interface.
    parser = argparse.ArgumentParser(prog='ARTDeco',
                                     description='Main script for Automatic ReadThrough DEteCtiOn-ARTDeco')
    parser.add_argument('-mode',
                        help='Mode in which to run ARTDeco. Options are preprocess, readthrough, diff_exp_read_in, \
                             get_dogs, and diff_exp_dogs.',action='store')
    parser.add_argument('-home-dir',help='Directory in which to run ARTDeco (default is current directory)',
                        action='store',default='.')
    parser.add_argument('-bam-files-dir',
                        help='Directory in which the BAM files are located (default is current directory)',
                        action='store',default='.')
    parser.add_argument('-layout',
                        help='Indicates whether the files are paired-end or single-end. Options are PE or SE.',
                        action='store')
    parser.add_argument('-stranded',
                        help='Indicates whether the files are stranded or unstranded. Options are True or False.',
                        action='store')
    parser.add_argument('-orientation',
                        help='Indicates whether the files are forward- or reverse-stranded. Options are Forward or '+
                             'Reverse. Required for stranded data',action='store')
    parser.add_argument('-single',
                        help='Indicates whether you want to create tag directories with a single file (useful for new'+\
                             'assemblies with lots of scaffolds).',default=False,action='store_true')
    parser.add_argument('-gtf-file',help='GTF file',action='store')
    parser.add_argument('-cpu',help='Maximum CPUs to use',action='store',type=int,default=1)
    parser.add_argument('-chrom-sizes-file',help='Chromosome sizes file')
    parser.add_argument('-read-in-dist',help='Read-in distance. Default is 1 kb.',type=int,default=1000)
    parser.add_argument('-readthrough-dist',help='Readthrough distance. Default is 10 kb.',type=int,default=10000)
    parser.add_argument('-intergenic-min-len',help='Minimum length for intergenic regions. Default is 100 bp.',type=int,
                        default=100)
    parser.add_argument('-intergenic-max-len',help='Maximum length for intergenic regions. Default is 15 kb.',type=int,
                        default=15000)
    parser.add_argument('-read-in-threshold', help='Threshold for considering read-in gene. Default is -1.',type=float,
                        default=-1)
    parser.add_argument('-read-in-fpkm',help='Minimum FPKM value for considering a gene. Default is 0.25 FPKM.',
                        type=float,default=0.25)
    parser.add_argument('-summary-genes', help='Number of genes for summarizing readthrough levels. Default is 1000.',
                        type=int,default=1000)
    parser.add_argument('-overwrite',help='Indicate whether to overwrite existing files',default=False,
                        action='store_true')
    parser.add_argument('-meta-file',help='Meta file',action='store',type=str)
    parser.add_argument('-comparisons-file',help='Comparisons file',type=str)
    parser.add_argument('-log2FC', help='Minimum log2 fold change for considering a gene upregulated. Default is 2.',
                        type=float,default=2)
    parser.add_argument('-pval', help='Maximum p-value for considering a gene upregulated. Default is 0.05.',type=float,
                        default=0.05)
    parser.add_argument('-min-dog-len',help='Minimum DoG length. Default is 4 kb.',type=int,default=4000)
    parser.add_argument('-dog-window',help='DoG window size. Default is 500 bp.',type=int,default=500)
    parser.add_argument('-min-dog-coverage',help='Minimum FPKM for DoG discovery. Default is 0.15 FPKM.',type=float,
                        default=0.15)
    parser.add_argument('-gene-types',help='Limit gene sets for reporting. Default is all gene types.',nargs='+',
                        action="store")
    parser.add_argument('-skip-bam-summary',help='Skip summary of BAM files (useful for time saving).',default=False,
                        action='store_true')
    args = parser.parse_known_args()[0]

    if args.mode in ['preprocess','readthrough','get_dogs','diff_exp_read_in','diff_exp_dogs']:
        print(f'Running {args.mode} mode...')
    else:
        print('No valid run mode specified... Will generate all files...')
        args.mode = None

    #Check if home and BAM file directories exist. If they do, load ARTDeco file structure.
    if os.path.isdir(args.home_dir) and os.path.isdir(args.bam_files_dir):
        print('Loading ARTDeco file structure...')
        artdeco_dir = ARTDecoDir(args.bam_files_dir,args.home_dir)
    elif not os.path.isdir(args.home_dir):
        print('User-specified home directory does not exist... Exiting...')
        sys.exit(1)
    else:
        print('User-specified BAM file directory does not exist... Exiting...')
        sys.exit(1)

    #Check for BAM files. If there are no BAM files, exit.
    if len(artdeco_dir.bam_files) == 0:
        print('No BAM files... Exiting...')
        sys.exit(1)

    #Overwrite specified.
    if args.overwrite:
        print('Overwrite specified... Will regenerate all files...')

    #Create summary file directory if it does not exist.
    if not os.path.isdir(artdeco_dir.summary_dir):
        os.mkdir(artdeco_dir.summary_dir)

    #Create preprocess_files directory if it does not exist.
    if not os.path.isdir(artdeco_dir.preprocess_dir):
        os.mkdir(artdeco_dir.preprocess_dir)

    #Generate meta and comparisons if it is needed.
    if (args.mode and (args.mode in ['diff_exp_read_in','diff_exp_dogs']) or
        (args.mode == 'preprocess' and args.meta_file)) or (not args.mode and args.meta_file):
        from .DESeq2 import reformat_meta,reformat_comparisons,generate_comparisons,load_deseq_dataset,run_deseq,\
            deseq_results

        #Specify whether meta file needs to be generated.
        if args.overwrite:
            new_meta = True
        elif os.path.isfile(artdeco_dir.meta_file):
            print('Reformatted meta file exists...')
            new_meta = False
        else:
            new_meta = True

        #Generate meta.
        if new_meta:

            #Check meta file.
            if args.meta_file and os.path.isfile(args.meta_file):

                meta = open(args.meta_file).readlines()
                meta_format = True
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
                    i += 1

                if meta_format:
                    print('Meta file properly formatted... Generating reformatted meta...')
                    reformat_meta(args.meta_file,artdeco_dir.preprocess_dir)
                else:
                    print('Meta file not properly formatted... Exiting...')
                    sys.exit(1)

            elif args.meta_file:
                print('Meta file does not exist... Exiting...')
                sys.exit(1)
            else:
                print('No meta file supplied... Exiting...')
                sys.exit(1)

        #Specify whether comparisons file needs to be generated.
        if args.overwrite:
            new_comparisons = True
        elif new_meta:
            new_comparisons = True
        elif os.path.isfile(artdeco_dir.comparisons_file):
            print('Reformatted comparisons file exists...')
            new_comparisons = False
        else:
            new_comparisons = True

        #Generate comparisons.
        if new_comparisons:

            #Grab groups.
            groups = [line.strip().split('\t')[1] for line in open(artdeco_dir.meta_file).readlines()[1:]]

            #Check comparisons file.
            if args.comparisons_file and os.path.isfile(args.comparisons_file):

                print('Comparisons file exists...')

                #Check format.
                comparisons = [line.strip().split('\t') for line in open(args.comparisons_file).readlines()]
                comparisons_lens = [len(line) for line in comparisons]

                #Check if lines are tab-separated formatted.
                if len(set(comparisons_lens)) == 1 and len(comparisons[0]) == 2:

                    comparisons_format = True
                    for line in comparisons:
                        line[0] = line[0].replace('-','_').replace(' ','_')
                        line[1] = line[1].replace('-', '_').replace(' ', '_')
                        if not line[0] in groups or not line[1] in groups:
                            comparisons_format = False
                else:
                    comparisons_format = False

                #If the file is properly formatted, reformat it. Otherwise, generate an all-by-all file.
                if comparisons_format:
                    print('Comparisons file properly formatted... Generating reformatted comparisons...')
                    reformat_comparisons(args.comparisons_file,artdeco_dir.preprocess_dir)
                else:
                    print('Comparisons file not properly formatted... Generating all-by-all comparisons file...')
                    generate_comparisons(artdeco_dir.meta_file,artdeco_dir.preprocess_dir)

            else:
                print('Comparison file does not exist or not provided... Generating comparisons file...')
                generate_comparisons(artdeco_dir.meta_file,artdeco_dir.preprocess_dir)

        comparisons = [line.strip().split('\t') for line in open(artdeco_dir.comparisons_file).readlines()]

    #Update files in ARTDeco directory.
    if os.path.exists(artdeco_dir.meta_file):
        artdeco_dir.set_diff_exp_output()

    #Generate output files.
    out_files = artdeco_dir.get_files(args.mode,os.path.exists(artdeco_dir.meta_file),args.overwrite)
    if out_files:
        print('ARTDeco will generate the following files:\n'+'\n'.join(out_files))
    else:
        print('All necessary files generated... Exiting...')
        sys.exit(1)

    #Update file structure.
    artdeco_dir.update_dir_lists(out_files)

    #Check if GTF file is needed for generating output.
    if artdeco_dir.gtf_needed:
        print('GTF file needed... Checking...')
        if args.gtf_file and os.path.isfile(args.gtf_file):
            print('GTF file exists...')
        elif args.gtf_file:
            print('User-supplied GTF file does not exist... Exiting...')
            sys.exit(1)
        else:
            print('No GTF file supplied... Exiting...')
            sys.exit(1)

    #Check if BAM file formats are needed. If they are, check if the user has specified them. If the user has not,
    #infer those formats. Summarize file if the format is needed.
    if artdeco_dir.format_needed:

        print('BAM file format needed... Checking... Will infer if not user-specified.')

        if args.layout:
            if str.lower(args.layout) in ['pe','se']:
                infer_layout = False
                if str.lower(args.layout) == 'pe':
                    print('BAM files specified as paired-end...')
                    pe = True
                else:
                    print('BAM files specified as single-end...')
                    pe = False
            else:
                print('Improper layout specified... Will infer...')
                infer_layout = True
        else:
            print('No layout specified... Will infer...')
            infer_layout = True

        if args.stranded:
            if str.lower(args.stranded) in ['true','false']:
                infer_strandedness = False
                if str.lower(args.stranded) == 'true':
                    print('BAM files specified as stranded...')
                    stranded = True
                else:
                    print('BAM files specified as unstranded...')
                    stranded = False
            else:
                print('Improper indication of strandedness...')
                infer_strandedness = True
        elif args.orientation and str.lower(args.orientation) in ['forward','reverse']:
            print('No strandedness specified but strand orientation specified... Will assign data as stranded...')
            infer_strandedness = False
            stranded = True
        else:
            print('No strandedness specified... Will infer...')
            infer_strandedness = True

        if args.orientation:
            if str.lower(args.orientation) in ['forward','reverse']:
                infer_orientation = False
                if str.lower(args.orientation) == 'forward':
                    print('BAM files specified as forward-strand oriented...')
                    flip = False
                else:
                    print('BAM files specified as reverse-strand oriented...')
                    flip = True
            else:
                print('Improper strand orientation specified... Will infer...')
                infer_orientation = True
        elif not infer_strandedness:
            if stranded:
                print('No strand orientation specified... Data is stranded... Will infer orientation...')
                infer_orientation = True
            else:
                print('No strand orientation specified... Data is unstranded... No need to infer orientation...')
                infer_orientation = False
                flip = False
        else:
            print('No strand orientation specified... Will infer...')
            infer_orientation = True

        #Infer layout if necessary.
        if infer_layout or infer_strandedness or infer_orientation:
            print('Will infer BAM formats...')

            #Check if full genes BED exists. If it doesn't, regenerate it.
            if os.path.isfile(artdeco_dir.genes_full) and not args.overwrite:
                print('Full genes BED file exists...')
            else:
                print('Generating full genes BED file...')
                parse_gtf(args.gtf_file,args.home_dir)
                out_files -= {artdeco_dir.genes_full,artdeco_dir.genes_condensed,artdeco_dir.gene_to_transcript,
                              artdeco_dir.gene_types}
                artdeco_dir.update_dir_lists(out_files)

            #Check file formats.
            print('Inferring BAM file formats...')
            formats = infer_experiments_group(artdeco_dir.bam_files,artdeco_dir.genes_full,
                                              min(args.cpu,len(artdeco_dir.bam_files)))

            #Check layout.
            if infer_layout:
                if len(set(x[1] for x in formats)) == 1:
                    pe = formats[0][1]
                    if pe:
                        print('All BAM files inferred as Paired-End...')
                    else:
                        print('All BAM files inferred as Single-End...')
                else:
                    print('Error... One or more files do not match in inferred format... Exiting...')
                    for f in formats:
                        out_str = f'BAM file {f[0]} inferred as '+output_inferred_format(f)
                        print(out_str)
                    sys.exit(1)

            #Check strandedness.
            if infer_strandedness:
                if len(set(x[2] for x in formats)) == 1:
                    stranded = formats[0][2]
                    if stranded:
                        print('All BAM files inferred as strand-specific...')
                    else:
                        print('All BAM files inferred as single-stranded...')
                else:
                    print('Error... One or more files do not match in inferred format... Exiting...')
                    for f in formats:
                        out_str = f'BAM file {f[0]} inferred as '+output_inferred_format(f)
                        print(out_str)
                    sys.exit(1)

            #Check strand orientation.
            if infer_orientation:
                if len(set(x[3] for x in formats)) == 1:
                    flip = formats[0][3]
                    if flip:
                        print('All BAM files inferred as reverse-strand oriented...')
                    else:
                        print('All BAM files inferred as forward-strand oriented...')
                else:
                    print('Error... One or more files do not match in inferred format... Exiting...')
                    for f in formats:
                        out_str = f'BAM file {f[0]} inferred as '+output_inferred_format(f)
                        print(out_str)
                    sys.exit(1)

        #Summarize files.
        if args.skip_bam_summary:
            print('Skipping summary of BAM file stats...')
        else:
            print('Summarizing BAM file stats...')
            summary_file = os.path.join(artdeco_dir.summary_dir,'bam_summary.txt')
            if os.path.isfile(summary_file):
                os.remove(summary_file)
            summary = summarize_bam_files(artdeco_dir.bam_files,args.cpu,pe,stranded,flip)
            for line in summary.split('\n'):
                print(line)
            with open(summary_file,'w') as f:
                f.write(summary)

    #Load chromosome sizes file if needed.
    if artdeco_dir.chrom_sizes_needed:
        if args.chrom_sizes_file and os.path.isfile(args.chrom_sizes_file):
            chrom_sizes_file = open(args.chrom_sizes_file)
            chrom_sizes = {}
            for line in chrom_sizes_file.readlines():
                line = line.strip().split('\t')
                chrom_sizes[line[0]] = int(line[1])
        elif args.chrom_sizes_file:
            print('Chromosome sizes file does not exist... Exiting...')
            sys.exit(1)
        else:
            print('No chromosome sizes file supplied... Exiting...')
            sys.exit(1)

    #Generate preprocessing files.
    if artdeco_dir.preprocessing_files:

        if not os.path.isdir(artdeco_dir.preprocess_dir):
            os.mkdir(artdeco_dir.preprocess_dir)

        if artdeco_dir.preprocessing_files & {artdeco_dir.genes_condensed,artdeco_dir.genes_full,
                                              artdeco_dir.gene_to_transcript,artdeco_dir.gene_types}:
            parse_gtf(args.gtf_file,args.home_dir)

        #Load genes.
        genes = pd.read_csv(artdeco_dir.genes_condensed,sep='\t',header=None,
                                         names=['Chrom','Start','Stop','Name','Score','Strand'])
        del genes['Score']
        genes = genes[['Name','Chrom','Strand','Start','Stop']]

        #Create read-in BED file.
        if artdeco_dir.read_in_bed in artdeco_dir.preprocessing_files:
            print('Generating read-in region BED file...')
            if stranded:
                read_in_df = create_stranded_read_in_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                        min_len=args.intergenic_min_len,upstream_dist=args.read_in_dist)
            else:
                read_in_df = create_unstranded_read_in_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                          min_len=args.intergenic_min_len,
                                                          upstream_dist=args.read_in_dist)
            #Format and create BED file.
            read_in_df['Score'] = 0
            read_in_df = read_in_df[['Chrom','Start','Stop','Name','Score','Strand']]
            read_in_df['Start'] = read_in_df['Start'].astype(int)
            read_in_df['Stop'] = read_in_df['Stop'].astype(int)

            #Output BED file.
            read_in_df.to_csv(artdeco_dir.read_in_bed,sep='\t',header=False,index=False)

        #Create readthrough BED file.
        if artdeco_dir.readthrough_bed in artdeco_dir.preprocessing_files:
            print('Generating readthrough region BED file...')
            if stranded:
                readthrough_df = create_stranded_downstream_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                               min_len=args.intergenic_min_len,
                                                               downstream_dist=args.readthrough_dist)
            else:
                readthrough_df = create_unstranded_downstream_df(genes,chrom_sizes,max_len=args.intergenic_max_len,
                                                                 min_len=args.intergenic_min_len,
                                                                 downstream_dist=args.readthrough_dist)

            readthrough_df['Score'] = 0
            readthrough_df = readthrough_df[['Chrom','Start','Stop','Name','Score','Strand']]
            readthrough_df['Start'] = readthrough_df['Start'].astype(int)
            readthrough_df['Stop'] = readthrough_df['Stop'].astype(int)

            readthrough_df.to_csv(artdeco_dir.readthrough_bed,sep='\t',header=False,index=False)

    #Creating tag directories.
    if artdeco_dir.tag_dirs:
        print('Creating tag directories...')
        make_multi_tag_dirs([artdeco_dir.tag_dir_to_bam[tag_dir] for tag_dir in artdeco_dir.tag_dirs],
                            artdeco_dir.preprocess_dir,flip,pe,stranded,args.single,min(len(artdeco_dir.tag_dirs),args.cpu))

    #Generate quantification files.
    if artdeco_dir.quantification_files:

        if not os.path.isdir(artdeco_dir.quantification_dir):
            print('Creating quantification directory...')
            os.mkdir(artdeco_dir.quantification_dir)

        if artdeco_dir.quantification_files & {artdeco_dir.gene_fpkm,artdeco_dir.gene_raw}:
            print('Generating gene expression files...')
            get_multi_gene_exp(artdeco_dir.all_tag_dirs,args.gtf_file,stranded,artdeco_dir.quantification_dir,args.cpu)

        if artdeco_dir.max_isoform in artdeco_dir.quantification_files:
            print('Getting maximum isoform...')
            get_max_isoform(artdeco_dir.gene_fpkm,artdeco_dir.gene_to_transcript,artdeco_dir.quantification_dir)

        if artdeco_dir.read_in_exp in artdeco_dir.quantification_files:
            print('Generating read-in expression file...')
            get_regions_exp((artdeco_dir.all_tag_dirs,artdeco_dir.read_in_bed,stranded,'-raw',
                             artdeco_dir.quantification_dir,min(len(artdeco_dir.all_tag_dirs),args.cpu)))

        if artdeco_dir.readthrough_exp in artdeco_dir.quantification_files:
            print('Generating readthrough expression file...')
            get_regions_exp((artdeco_dir.all_tag_dirs, artdeco_dir.readthrough_bed,stranded,'-raw',
                             artdeco_dir.quantification_dir,min(len(artdeco_dir.all_tag_dirs),args.cpu)))

    #Generate readthrough files.
    if artdeco_dir.readthrough_files:

        if not os.path.isdir(artdeco_dir.readthrough_dir):
            print('Creating readthrough directory...')
            os.mkdir(artdeco_dir.readthrough_dir)

        if artdeco_dir.read_in_levels in artdeco_dir.readthrough_files:
            print('Generate read-in vs. expression file...')
            get_gene_v_intergenic(artdeco_dir.gene_raw,artdeco_dir.gene_fpkm,artdeco_dir.max_isoform,
                                  artdeco_dir.read_in_exp,'Read-In',artdeco_dir.read_in_levels)

        if artdeco_dir.corrected_exp in artdeco_dir.readthrough_files:
            print('Correcting gene expression using read-in information...')
            deconvolute_exp(artdeco_dir.read_in_levels,artdeco_dir.corrected_exp)

        if artdeco_dir.readthrough_levels in artdeco_dir.readthrough_files:
            print('Generate readthrough vs. expression file...')
            get_gene_v_intergenic(artdeco_dir.gene_raw,artdeco_dir.gene_fpkm,artdeco_dir.max_isoform,
                                  artdeco_dir.readthrough_exp,'Readthrough',artdeco_dir.readthrough_levels)

        if artdeco_dir.read_in_assignments in artdeco_dir.readthrough_files:
            print(f'Read-in genes assigned with read-in level threshold is {args.read_in_threshold} and read-in FPKM '+\
                  f'threshold is {args.read_in_fpkm}...')
            if args.gene_types and len(open(artdeco_dir.gene_types).readlines()) > 1:
                print('Using the following gene types: ' + ', '.join(args.gene_types) + '...')
            else:
                args.gene_types = None
                print('Using all genes...')
            assign_genes(artdeco_dir.read_in_levels,args.read_in_threshold,args.read_in_fpkm,
                         artdeco_dir.read_in_assignments,artdeco_dir.gene_types,args.gene_types)

        #Summarize output.
        print('Summarizing readthrough output...')

        if args.gene_types and len(open(artdeco_dir.gene_types).readlines()) > 1:
            print('Using the following gene types: '+', '.join(args.gene_types)+'...')
        else:
            args.gene_types = None
            print('Using all genes...')

        summary_file = os.path.join(artdeco_dir.summary_dir,'readthrough_summary.txt')

        if os.path.isfile(summary_file):
            os.remove(summary_file)

        expts = []
        for f in os.listdir(artdeco_dir.bam_files_dir):
            if f[-4:] == '.bam':
                expt = f[:-4].replace('-', '_').replace(' ', '_')
                expts.append(expt)

        summary = summarize_readthrough_stats(artdeco_dir.read_in_levels,expts,'Read-In',
                                              args.summary_genes,artdeco_dir.gene_types,args.gene_types)

        if os.path.isfile(artdeco_dir.readthrough_levels):
            summary += '\n'+summarize_readthrough_stats(artdeco_dir.readthrough_levels,expts,'Readthrough',
                                                        args.summary_genes,artdeco_dir.gene_types,args.gene_types)

        summary += '\n'+summarize_read_in_assignments(artdeco_dir.read_in_assignments,expts,args.read_in_threshold,
                                                      args.read_in_fpkm)

        for line in summary.split('\n'):
            print(line)

        with open(summary_file,'w') as f:
            f.write(summary)

    #Generate DoG files.
    if artdeco_dir.dogs_files:

        if not os.path.isdir(artdeco_dir.dogs_dir):
            print('Creating DoG output directory...')
            os.mkdir(artdeco_dir.dogs_dir)

        if artdeco_dir.dogs_beds:

            print(f'Finding DoGs...\nGet genes with potential DoGs with minimum length of {args.min_dog_len} bp, a '+
                  f'minimum coverage of {args.min_dog_coverage} FPKM, and screening window of {args.dog_window} bp...')
            screening_genes = get_dog_screening(artdeco_dir.genes_condensed,args.min_dog_len)

            print(f'Generate initial screening BED file for DoGs with minimum length {args.min_dog_len} bp and window'+
                  f' size {args.dog_window} bp...')
            generate_screening_bed(screening_genes,args.min_dog_len,args.dog_window,args.home_dir)

            print(f'Initial screening coverage for DoGs with minimum length of {args.min_dog_len} bp...')

            #Screen for coverage threshold.
            dogs_tag_dirs = [artdeco_dir.dogs_bed_to_tagdir[dogs_bed] for dogs_bed in artdeco_dir.dogs_beds]
            screening_coverage_dfs = get_multi_interval_coverage(dogs_tag_dirs,args.home_dir,
                                                                 [os.path.join(args.home_dir,'intervals.bed')],
                                                                 args.min_dog_coverage,stranded,
                                                                 min(args.cpu,len(dogs_tag_dirs)))

            #Find genes that pass threshold for minimum length.
            screening_genes_dfs = {}
            for i in range(len(screening_coverage_dfs)):
                new_df = screening_coverage_dfs[i].copy()
                new_df = new_df.groupby('Name').count()
                new_df = new_df[new_df.ID == args.min_dog_len/args.dog_window]
                screening_genes_dfs[dogs_tag_dirs[i]] = screening_genes[screening_genes.Name.isin(new_df.index)].copy()

            #Remove screening BED files.
            os.remove(os.path.join(args.home_dir,'intervals.bed'))

            print('Generate screening BED file for pre-screened DoGs...')

            generate_full_screening_bed(dogs_tag_dirs,artdeco_dir.genes_condensed,screening_genes_dfs,
                                        artdeco_dir.read_in_assignments,args.chrom_sizes_file,args.dog_window,args.cpu,
                                        args.home_dir)

            print('Screening coverage for pre-screened DoGs...')

            #Screen for coverage threshold.
            screening_coverage_dfs = get_multi_interval_coverage(dogs_tag_dirs,args.home_dir,
                                                                 [os.path.join(args.home_dir,
                                                                               tag_dir.split('/')[-1]+'.bed') for
                                                                  tag_dir in dogs_tag_dirs],args.min_dog_coverage,
                                                                 stranded,min(args.cpu,len(dogs_tag_dirs)))

            #Remove screening BED files.
            for tag_dir in dogs_tag_dirs:
                expt_name = tag_dir.split('/')[-1]
                os.remove(os.path.join(args.home_dir,f'{expt_name}.bed'))

            print('Discovering DoG coordinates for pre-screened DoGs and output BED files...')
            get_multi_dog_beds(screening_genes,screening_coverage_dfs,args.dog_window,args.cpu,dogs_tag_dirs,
                               artdeco_dir.dogs_dir)

        if artdeco_dir.all_dogs_bed in artdeco_dir.dogs_files or artdeco_dir.dogs_beds:
            print('Merge DoGs into a single annotation...')
            merge_dogs(artdeco_dir.all_dogs_beds,artdeco_dir.dogs_dir)

        if artdeco_dir.dogs_raw & artdeco_dir.dogs_files or artdeco_dir.dogs_fpkm & artdeco_dir.dogs_files:
            print('Generating expression data for DoGs in individual experiments...')
            get_dog_exp(artdeco_dir.all_tag_dirs,artdeco_dir.all_dogs_beds,stranded,artdeco_dir.dogs_dir,args.cpu)

        if artdeco_dir.all_dogs_fpkm in artdeco_dir.dogs_files or artdeco_dir.dogs_beds:
            print('Generating expression data in FPKM for all DoGs...')
            get_regions_exp((artdeco_dir.all_tag_dirs,artdeco_dir.all_dogs_bed,stranded,'-fpkm',artdeco_dir.dogs_dir,
                             min(args.cpu,len(artdeco_dir.all_tag_dirs))))

        if artdeco_dir.all_dogs_raw in artdeco_dir.dogs_files or artdeco_dir.dogs_beds:
            print('Generating raw expression data for all DoGs...')
            get_regions_exp((artdeco_dir.all_tag_dirs,artdeco_dir.all_dogs_bed,stranded,'-raw',artdeco_dir.dogs_dir,
                             min(args.cpu,len(artdeco_dir.all_tag_dirs))))

        #Summarize DoG files.
        summary_file = os.path.join(artdeco_dir.summary_dir,'dogs_summary.txt')
        if os.path.isfile(summary_file):
            os.remove(summary_file)

        summary = summarize_all_dogs(artdeco_dir.all_dogs_bed,artdeco_dir.all_dogs_beds,artdeco_dir.all_dogs_fpkm,
                                     artdeco_dir.all_dogs_fpkm_expts,args.min_dog_len,args.min_dog_coverage,
                                     args.dog_window)

        for line in summary.split('\n'):
            print(line)

        with open(summary_file,'w') as f:
            f.write(summary)

    #Generate differential expression output.
    try:
        if artdeco_dir.diff_exp_files:

            if not os.path.isdir(artdeco_dir.diff_exp_dir):
                print('Creating differential expression output directory...')
                os.mkdir(artdeco_dir.diff_exp_dir)

            print('Running DESeq2 on gene expression data...')
            dds = load_deseq_dataset(artdeco_dir.gene_raw,artdeco_dir.meta_file)
            dds_results = run_deseq(dds)

            #Output results.
            print('Output DESeq2 results...')
            for condition1, condition2 in comparisons:
                deseq_results(dds_results,condition1,condition2,artdeco_dir.diff_exp_dir)
    except:
        pass

    #Generate differential expression with read-in information.
    try:
        if artdeco_dir.diff_exp_read_in_files:

            if not os.path.isdir(artdeco_dir.diff_exp_read_in_dir):
                print('Creating differential expression with read-in information directory...')
                os.mkdir(artdeco_dir.diff_exp_read_in_dir)

            #Join differential expression information to read-in information. Use this to infer read-in genes.
            print('Combining differential expression results and read-in information... Inferring read-in genes for '+\
                  f'upregulated genes with log2 fold change > {args.log2FC}, p-value < {args.pval}, and FPKM > '+\
                  f'{args.read_in_fpkm}... Read-in level threshold is {args.read_in_threshold}...')
            if args.gene_types and len(open(artdeco_dir.gene_types).readlines()) > 1:
                print('Using the following gene types: ' + ', '.join(args.gene_types) + '...')
            else:
                args.gene_types = None
                print('Using all genes...')
            for condition1,condition2 in comparisons:

                read_in_diff_exp(artdeco_dir.read_in_levels,artdeco_dir.meta_file,
                                 os.path.join(artdeco_dir.diff_exp_dir,f'{condition1}-{condition2}-results.txt'),
                                 artdeco_dir.diff_exp_read_in_dir)

                assign_read_in_genes(os.path.join(artdeco_dir.diff_exp_read_in_dir,
                                                  f'{condition1}-{condition2}-read_in.txt'),args.log2FC,args.pval,
                                     args.read_in_fpkm,args.read_in_threshold,artdeco_dir.gene_types,args.gene_types,
                                     artdeco_dir.diff_exp_read_in_dir)

            #Summarize output.
            print('Summarize read-in gene inference with differential expression information...')
            summary_file = os.path.join(artdeco_dir.summary_dir,'diff_exp_read_in_summary.txt')
            if os.path.isfile(summary_file):
                os.remove(summary_file)

            assignment_files = []
            for condition1,condition2 in comparisons:
                assignment_files.append(os.path.join(artdeco_dir.diff_exp_read_in_dir,
                                                     f'{condition1}-{condition2}-read_in_assignment.txt'))

            summary = summarize_diff_exp_read_in_assignments(assignment_files,args.log2FC,args.pval,args.read_in_fpkm,
                                                             args.read_in_threshold)

            for line in summary.split('\n'):
                print(line)

            with open(summary_file,'w') as f:
                f.write(summary)

    except:
        pass

    #Generate differential expression output for DoGs.
    try:
        if artdeco_dir.diff_exp_dogs_files:

            if not os.path.isdir(artdeco_dir.diff_exp_dogs_dir):
                print('Creating differential expression for DoGs directory...')
                os.mkdir(artdeco_dir.diff_exp_dogs_dir)

            print('Running DESeq2 on DoGs...')
            dds = load_deseq_dataset(artdeco_dir.all_dogs_raw,artdeco_dir.meta_file)
            dds_results = run_deseq(dds)

            #Output results.
            print('Output DESeq2 results...')
            for condition1,condition2 in comparisons:
                deseq_results(dds_results,condition1,condition2,artdeco_dir.diff_exp_dogs_dir)
    except:
        pass

if __name__ == '__main__':
    main()