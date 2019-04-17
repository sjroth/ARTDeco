'''
Script that can run the intergenic mode.
'''
from ARTDeco.modules.misc import infer_experiments_group,output_inferred_format,get_regions_exp
from ARTDeco.modules.intergenic import get_multi_gene_exp,get_max_isoform,get_gene_v_intergenic,assign_genes

import sys
import os

def main(argv):

    home_dir = argv[0]
    gtf_file = argv[1]
    cpu = int(argv[2])
    read_in_threshold = float(argv[3])
    read_in_fpkm = float(argv[4])
    overwrite = bool(argv[5] == 'True')

    #Grab tag directories based upon input bam files
    tag_dirs = []
    bam_files = []
    for f in os.listdir(home_dir):
        if f[-4:] == '.bam':
            tag_dirs.append(os.path.join(home_dir,'preprocess_files',f[:-4]))
            bam_files.append(f)

    #Check that the prerequisites exist.
    print('Checking prerequisites...')
    if os.path.isdir(os.path.join(home_dir,'preprocess_files')):
        if os.path.isfile(os.path.join(home_dir,'preprocess_files','genes.full.bed')) and \
                os.path.isfile(os.path.join(home_dir,'preprocess_files','gene_to_transcript.txt')):
            print('Necessary gene annotation files exist...')
        else:
            print('Necessary gene annotation files do not exist... Please run preprocessing mode... Exiting...')
            sys.exit(1)

        if os.path.isfile(os.path.join(home_dir,'preprocess_files','read_in.bed')) and \
                os.path.isfile(os.path.join(home_dir,'preprocess_files','readthrough.bed')):
            print('Intergenic BED files exist...')
        else:
            print('At least one of the intergenic BED files does not exist... Please run preprocessing mode... \
            Exiting...')

        for tag_dir in tag_dirs:

            if not os.path.isdir(tag_dir):
                print('Tag directories missing... Please run preprocessing mode... Exiting...')
                sys.exit(1)
        print('All tag directories exist...')

    else:
        print('Preprocessing files do not exist... Please run preprocessing mode... Exiting...')
        sys.exit(1)

    #Check if user-supplied GTF file exists.
    if os.path.isfile(gtf_file):
        print('GTF file exists...')
    else:
        print('Invalid GTF file supplied... Exiting...')
        sys.exit(1)

    #Create quantification directory.
    if os.path.isdir(os.path.join(home_dir,'quantification')):
        print('Quantification directory exists...')
    else:
        print('Creating quantification directory...')
        os.mkdir(os.path.join(home_dir,'quantification'))

    #Create quantification directory.
    if os.path.isdir(os.path.join(home_dir,'intergenic')):
        print('Intergenic vs. gene expression directory exists...')
    else:
        print('Creating intergenic vs. gene expression directory...')
        os.mkdir(os.path.join(home_dir,'intergenic'))

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

    #Check if gene expression and maximum isoform files exist. If one of them does not, create all three.
    if os.path.isfile(os.path.join(home_dir,'quantification','gene.exp.raw.txt')) and \
            os.path.isfile(os.path.join(home_dir,'quantification','gene.exp.fpkm.txt')) and \
            os.path.isfile(os.path.join(home_dir,'quantification','max_isoform.txt')) and not overwrite:
        print('Gene expression files exist...')
    else:
        print('Generating gene expression and maximum isoform files...')
        get_multi_gene_exp(tag_dirs,gtf_file,stranded,os.path.join(home_dir,'quantification'),cpu)
        get_max_isoform(os.path.join(home_dir,'quantification','gene.exp.fpkm.txt'),
                        os.path.join(home_dir,'preprocess_files','gene_to_transcript.txt'),
                        os.path.join(home_dir,'quantification'))

    #Check if read-in and readthrough expression files exist. If one of them does not, create both of them.
    if os.path.isfile(os.path.join(home_dir,'quantification','read_in.exp.txt')) and \
            os.path.isfile(os.path.join(home_dir,'quantification','readthrough.exp.txt')) and not overwrite:
        print('Intergenic expression files exist...')
    else:
        print('Generating intergenic expression files...')
        get_regions_exp((tag_dirs,os.path.join(home_dir,'preprocess_files','read_in.bed'),stranded,'-raw',
                         os.path.join(home_dir,'quantification'),min(len(tag_dirs),cpu)))
        get_regions_exp((tag_dirs,os.path.join(home_dir,'preprocess_files','readthrough.bed'),stranded,'-raw',
                         os.path.join(home_dir,'quantification'),min(len(tag_dirs),cpu)))

    #Check if intergenic vs. gene expression files exist. IF one of them does not, create both of them.
    if os.path.isfile(os.path.join(home_dir,'intergenic','read_in.txt')) and \
            os.path.isfile(os.path.join(home_dir,'intergenic','readthrough.txt')) and \
            os.path.isfile(os.path.join(home_dir,'intergenic','read_in_assignments.txt')) and not overwrite:
        print('Intergenic vs. gene expression files and read-in gene assignments exist...')
    else:
        print(f'Generating intergenic vs. expression files and read-in gene assignments... Read-in level threshold is'+\
              f' {read_in_threshold} and read-in FPKM threshold is {read_in_fpkm}...')

        get_gene_v_intergenic(os.path.join(home_dir,'quantification','gene.exp.raw.txt'),
                              os.path.join(home_dir,'quantification','gene.exp.fpkm.txt'),
                              os.path.join(home_dir,'quantification','max_isoform.txt'),
                              os.path.join(home_dir,'quantification','read_in.raw.txt'),'Read-In',
                              os.path.join(home_dir,'intergenic','read_in.txt'))

        get_gene_v_intergenic(os.path.join(home_dir,'quantification','gene.exp.raw.txt'),
                              os.path.join(home_dir,'quantification','gene.exp.fpkm.txt'),
                              os.path.join(home_dir,'quantification','max_isoform.txt'),
                              os.path.join(home_dir,'quantification','readthrough.raw.txt'),'Readthrough',
                              os.path.join(home_dir,'intergenic','readthrough.txt'))

        assign_genes(os.path.join(home_dir,'intergenic','read_in.txt'),read_in_threshold,read_in_fpkm,
                     os.path.join(home_dir,'intergenic','read_in_assignments.txt'))

if __name__ == '__main__':
    main(sys.argv[1:])