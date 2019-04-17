import argparse
import os
import sys
import subprocess

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
    parser.add_argument('-read-in-dist',help='Read-in distance',type=int,default=1000)
    parser.add_argument('-readthrough-dist',help='Readthrough distance',type=int,default=5000)
    parser.add_argument('-intergenic-min-len',help='Minimum length for intergenic regions',type=int,default=100)
    parser.add_argument('-intergenic-max-len',help='Maximum length for intergenic regions',type=int,default=15000)
    parser.add_argument('-read-in-threshold', help='Threshold for considering read-in gene',type=float,default=0)
    parser.add_argument('-fpkm',help='Minimum FPKM value for considering a gene',type=float,default=0.25)
    parser.add_argument('-overwrite',help='Indicate whether to overwrite existing files',default=False,
                        action='store_true')
    parser.add_argument('-meta-file',help='Meta file',action='store',type=str,default='.')
    parser.add_argument('-comparisons-file',help='Comparisons file',type=str,default='.')
    parser.add_argument('-log2FC', help='Minimum log2 fold change for considering a gene upregulated',type=float,
                        default=2)
    parser.add_argument('-pval', help='Maximum p-value for considering a gene upregulated',type=float,default=0.05)
    parser.add_argument('-min-dog-len',help='Minimum DoG length',type=int,default=4000)
    parser.add_argument('-dog-window',help='DoG window size',type=int,default=500)
    parser.add_argument('-min_dog_coverage',help='Minimum FPKM for DoG discovery',type=float,default=0.1)
    args = parser.parse_known_args()[0]

    #Check home directory for BAM files.
    no_bams = True
    for f in os.listdir(args.home_dir):
        if f[-4:] == '.bam':
            no_bams = False
            break

    #If there are no BAM files, exit.
    if no_bams:
        print('No BAM files... Exiting...')
        sys.exit(1)

    #Preprocess mode.
    if args.mode == 'preprocess':
        subprocess.call(['python','preprocess.py',args.home_dir,args.gtf_file,str(args.cpu),args.chrom_sizes_file,
                         str(args.intergenic_max_len),str(args.intergenic_min_len),str(args.read_in_dist),
                         str(args.readthrough_dist),str(args.overwrite)])

    if args.mode == 'intergenic':
        subprocess.call(['python','intergenic.py',args.home_dir,args.gtf_file,str(args.cpu),str(args.read_in_threshold),
                         str(args.fpkm),str(args.overwrite)])

    if args.mode == 'diff_exp_read_in':
        subprocess.call(['python','diff_exp_read_in.py',args.home_dir,args.meta_file,args.comparisons_file,
                         str(args.log2FC),str(args.pval),str(args.fpkm),str(args.read_in_threshold),
                         str(args.overwrite)])

    if args.mode == 'get_dogs':
        subprocess.call(['python','get_dogs.py',args.home_dir,str(args.cpu),str(args.min_dog_len),str(args.dog_window),
                         str(args.min_dog_coverage),args.chrom_sizes_file,str(args.overwrite)])

    if args.mode == 'diff_exp_dogs':
        subprocess.call(['python','diff_exp_dogs.py',args.home_dir,args.meta_file,args.comparisons_file,
                         str(args.overwrite)])

if __name__ == '__main__':
    main()