'''
Script that can run the preprocessing step.
'''
from ARTDeco.modules.preprocess import parse_gtf,create_stranded_downstream_df,create_stranded_read_in_df,\
    create_unstranded_downstream_df,create_unstranded_read_in_df,make_multi_tag_dirs
from ARTDeco.modules.misc import infer_experiments_group,output_inferred_format

import sys
import os
import pandas as pd

def main(argv):

    home_dir = argv[0]
    gtf_file = argv[1]
    cpu = int(argv[2])
    chrom_sizes_file = argv[3]
    intergenic_max_len = int(argv[4])
    intergenic_min_len = int(argv[5])
    read_in_dist = int(argv[6])
    readthrough_dist = int(argv[7])
    overwrite = bool(argv[8] == 'True')

    #Check home directory for BAM files.
    bam_files = []
    for f in os.listdir(home_dir):
        if f[-4:] == '.bam':
            bam_files.append(f)

    #Check for existence of a gene annotation directory with all necessary files. Create file if necessary.
    if not os.path.isfile(os.path.join(home_dir,'preprocess_files','genes_condensed.bed')) or \
            not os.path.isfile(os.path.join(home_dir,'preprocess_files','gene_to_transcript.txt')) or overwrite:
        print('Incomplete gene annotation files or overwrite specified... Will create gene annotations...')

        #Create directory if it does not already exist.
        if not os.path.isdir(os.path.join(home_dir,'preprocess_files')):
            os.makedirs(os.path.join(home_dir,'preprocess_files'))

        #Check if user-supplied GTF file exists. If it does, create the necessary files.
        if os.path.isfile(gtf_file):
            parse_gtf(gtf_file=gtf_file,home_dir=home_dir)
        else:
            print('Invalid GTF file supplied... Exiting...')
            sys.exit(1)

    else:
        print('Gene annotation files exist...')

    #Infer BAM file formats.
    print('Inferring file formats...')
    formats = infer_experiments_group([os.path.join(home_dir,bam_file) for bam_file in bam_files],
                                      os.path.join(home_dir,'preprocess_files','genes.full.bed'),
                                      min(cpu,len(bam_files)))

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

    #Create read-in and downstream BED files.
    if not os.path.isfile(os.path.join(home_dir,'preprocess_files','read_in.bed')) or \
            not os.path.isfile(os.path.join(home_dir,'preprocess_files','readthrough.bed')) or overwrite:

        print('Creating intergenic BED files...')

        #Load genes.
        genes = pd.read_csv(os.path.join(home_dir,'preprocess_files/genes_condensed.bed'),sep='\t',header=None,
                            names=['Chrom','Start','Stop','Name','Score','Strand'])
        del genes['Score']
        genes = genes[['Name','Chrom','Strand','Start','Stop']]

        #Load chromosome sizes file.
        if os.path.isfile(chrom_sizes_file):
            chrom_sizes_file = open(chrom_sizes_file)
            chrom_sizes = {}
            for line in chrom_sizes_file.readlines():
                line = line.strip().split('\t')
                chrom_sizes[line[0]] = int(line[1])
        else:
            print('Chromosome sizes file does not exist... Exiting...')
            sys.exit(1)

        #Create intergenic dataframes.
        if stranded:
            read_in_df = create_stranded_read_in_df(genes,chrom_sizes,max_len=intergenic_max_len,
                                                    min_len=intergenic_min_len,upstream_dist=read_in_dist)
            readthrough_df = create_stranded_downstream_df(genes,chrom_sizes,max_len=intergenic_max_len,
                                                           min_len=intergenic_min_len,downstream_dist=readthrough_dist)
        else:
            read_in_df = create_unstranded_read_in_df(genes,chrom_sizes,max_len=intergenic_max_len,
                                                      min_len=intergenic_min_len,upstream_dist=read_in_dist)
            readthrough_df = create_unstranded_downstream_df(genes,chrom_sizes,max_len=intergenic_max_len,
                                                             min_len=intergenic_min_len,
                                                             downstream_dist=readthrough_dist)

        #Format and create BED files.
        read_in_df['Score'] = 0
        read_in_df = read_in_df[['Chrom','Start','Stop','Name','Score','Strand']]
        read_in_df['Start'] = read_in_df['Start'].astype(int)
        read_in_df['Stop'] = read_in_df['Stop'].astype(int)

        readthrough_df['Score'] = 0
        readthrough_df = readthrough_df[['Chrom','Start','Stop','Name','Score','Strand']]
        readthrough_df['Start'] = readthrough_df['Start'].astype(int)
        readthrough_df['Stop'] = readthrough_df['Stop'].astype(int)

        read_in_df.to_csv(os.path.join(home_dir,'preprocess_files','read_in.bed'),sep='\t',header=False,index=False)
        readthrough_df.to_csv(os.path.join(home_dir,'preprocess_files','readthrough.bed'),sep='\t',header=False,
                              index=False)

    else:
        print('Intergenic files already exist...')

    #Create tag directories if they don't exist.
    tag_dir_bams = []
    for bam_file in bam_files:
        tag_dir = os.path.join(home_dir,'preprocess_files',bam_file[:-4])
        if not os.path.isdir(tag_dir) or overwrite:
            tag_dir_bams.append(os.path.join(home_dir,bam_file))

    if tag_dir_bams:
        print('Creating tag directories...')
        make_multi_tag_dirs(tag_dir_bams,os.path.join(home_dir,'preprocess_files'),flip,pe,stranded,
                            min(len(bam_files),cpu))
    else:
        print('Tag directories already exist...')

if __name__ == '__main__':
    main(sys.argv[1:])