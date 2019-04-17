'''
Script that can run the intergenic mode.
'''
from ARTDeco.modules.misc import infer_experiments_group,output_inferred_format,get_regions_exp
from ARTDeco.modules.get_dogs import get_dog_screening,generate_screening_bed,get_multi_interval_coverage,\
    generate_full_screening_bed,get_multi_dog_beds,merge_dogs,get_dog_exp

import os
import sys

def main(argv):

    home_dir = argv[0]
    cpu = int(argv[1])
    min_dog_len = int(argv[2])
    dog_window = int(argv[3])
    min_dog_coverage = float(argv[4])
    chrom_sizes_file = argv[5]
    overwrite = bool(argv[6] == 'True')

    #Grab tag directories based upon input bam files
    tag_dirs = []
    bam_files = []
    for f in os.listdir(home_dir):
        if f[-4:] == '.bam':
            tag_dirs.append(os.path.join(home_dir,'preprocess_files',f[:-4]))
            bam_files.append(f)

    #Check that the prerequisites exist.
    print('Checking prerequisites...')
    if os.path.isfile(os.path.join(home_dir,'preprocess_files','genes_condensed.bed')) and \
            os.path.isfile(os.path.join(home_dir,'preprocess_files','genes.full.bed')):
                print('Necessary gene annotation files exist...')
    else:
            print('Necessary gene annotation files do not exist... Please run preprocessing mode... Exiting...')
            sys.exit(1)

    for tag_dir in tag_dirs:
        if not os.path.isdir(tag_dir):
            print('Tag directories missing... Please run preprocessing mode... Exiting...')
            sys.exit(1)
    print('All tag directories exist...')

    if os.path.isfile(os.path.join(home_dir,'intergenic','read_in_assignments.txt')):
        print('Read-in assignments exist...')
    else:
        print('Read-in assignments do not exist... Please run intergenic mode... Exiting...')
        sys.exit(1)

    #Infer formats of BAM files.
    print('Inferring file formats...')
    formats = infer_experiments_group([os.path.join(home_dir,bam_file) for bam_file in bam_files],
                                      os.path.join(home_dir,'preprocess_files','genes.full.bed'),
                                      min(cpu,len(bam_files)))

    #Check that all of the BAM files are of the same format.
    if len(set(x[1] for x in formats)) == 1 and len(set(x[2] for x in formats)) == 1 and \
            len(set(x[3] for x in formats)) == 1:

        pe = formats[0][1]
        stranded = formats[0][2]
        flip = formats[0][3]

        out_str = 'All BAM files are '+output_inferred_format(formats[0])

        print(out_str)

    else:

        print('Error... One or more files do not match in inferred format... Exiting...')

        for f in formats:
            out_str = f'BAM file {f[0]} inferred as '+output_inferred_format(f)
            print(out_str)

        sys.exit(1)

    #Create DoG directory if one does not exist.
    if os.path.isdir(os.path.join(home_dir,'dogs')):
        print('DoG directory exists...')
    else:
        print('Creating DoG directory...')
        os.mkdir(os.path.join(home_dir,'dogs'))

    dog_files = []
    all_bed_exist = True
    all_exp_exist = True
    for tag_dir in tag_dirs:

        dog = os.path.join(home_dir,'dogs',tag_dir.split('/')[-1])+'.dogs.'

        dog_files.append(dog+'bed')

        if not os.path.isfile(dog+'bed'):
            all_bed_exist = False

        if not os.path.isfile(dog+'raw.txt') or not os.path.isfile(dog+'fpkm.txt'):
            all_exp_exist = False

    if not os.path.isfile(os.path.join(home_dir,'dogs','all_dogs.bed')):
        all_bed_exist = False

    if not os.path.isfile(os.path.join(home_dir,'dogs','all_dogs.raw.txt')) or \
            not os.path.isfile(os.path.join(home_dir,'dogs','all_dogs.fpkm.txt')):
        all_exp_exist = False

    if not all_bed_exist or overwrite:
        print('Finding DoGs...')

        print(f'Get genes with potential DoGs with minimum length of {min_dog_len} bp...')
        screening_genes = get_dog_screening(os.path.join(home_dir,'preprocess_files','genes_condensed.bed'),min_dog_len)

        print(f'Generate initial screening BED file for DoGs with minimum length {min_dog_len} bp and window size '+\
              f'{dog_window} bp...')
        generate_screening_bed(screening_genes,min_dog_len,dog_window,home_dir)

        print(f'Initial screening coverage for DoGs with minimum length of {min_dog_len} bp...')

        #Screen for coverage threshold.
        screening_coverage_dfs = get_multi_interval_coverage(tag_dirs,home_dir,[os.path.join(home_dir,'intervals.bed')],
                                                             min_dog_coverage,stranded,min(cpu,len(tag_dirs)))

        #Find genes that pass threshold for minimum length.
        screening_genes_dfs = {}
        for i in range(len(screening_coverage_dfs)):
            new_df = screening_coverage_dfs[i].copy()
            new_df = new_df.groupby('Name').count()
            new_df = new_df[new_df.ID == min_dog_len/dog_window]
            screening_genes_dfs[tag_dirs[i]] = screening_genes[screening_genes.Name.isin(new_df.index)].copy()

        #Remove screening BED files.
        os.remove(os.path.join(home_dir,'intervals.bed'))

        print('Generate screening BED file for pre-screened DoGs...')

        #Check if chromosome sizes file exists.
        if not os.path.isfile(chrom_sizes_file):
            print('Chromosome sizes file does not exist... Exiting...')
            sys.exit(1)

        generate_full_screening_bed(tag_dirs,os.path.join(home_dir,'preprocess_files','genes_condensed.bed'),
                                    screening_genes_dfs,os.path.join(home_dir,'intergenic','read_in_assignments.txt'),
                                    chrom_sizes_file,dog_window,cpu,home_dir)

        print('Screening coverage for pre-screened DoGs...')

        #Screen for coverage threshold.
        screening_coverage_dfs = get_multi_interval_coverage(tag_dirs,home_dir,
                                        [os.path.join(home_dir,tag_dir.split('/')[-1]+'.bed') for tag_dir in tag_dirs],
                                                             min_dog_coverage,stranded,min(cpu,len(tag_dirs)))

        #Remove screening BED files.
        for tag_dir in tag_dirs:
            expt_name = tag_dir.split('/')[-1]
            os.remove(os.path.join(home_dir,f'{expt_name}.bed'))

        print('Discovering DoG coordinates for pre-screened DoGs and output BED files...')
        get_multi_dog_beds(screening_genes,screening_coverage_dfs,dog_window,cpu,tag_dirs,os.path.join(home_dir,'dogs'))

        print('Merge DoGs into a single annotation...')
        merge_dogs(dog_files,os.path.join(home_dir,'dogs'))

    else:
        print('All DoG annotations exist...')

    if not all_exp_exist or overwrite:

        print('Quantifying DoGs...')

        print('Quantify DoGs for each experiment...')
        get_dog_exp(tag_dirs,dog_files,stranded,os.path.join(home_dir,'dogs'),cpu)

        print('Quantify all DoGs for all experiments...')
        get_regions_exp((tag_dirs,os.path.join(home_dir,'dogs','all_dogs.bed'),stranded,'-raw',
                         os.path.join(home_dir,'dogs'),min(cpu,len(tag_dirs))))
        get_regions_exp((tag_dirs,os.path.join(home_dir,'dogs','all_dogs.bed'),stranded,'-fpkm',
                         os.path.join(home_dir,'dogs'),min(cpu,len(tag_dirs))))
    else:
        print('DoG quantification files exist...')

if __name__ == '__main__':
    main(sys.argv[1:])