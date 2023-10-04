def split_by_chrom_helper(vcf_filename):
    '''helper to split vcf file rows by chromosomes'''
    out_dict = defaultdict(list)
    with open(vcf_filename,'r') as vcf_file:
        for line in vcf_file.readlines():
            vcf_line = line.rstrip().split('\t')
            if '_' not in vcf_line[0] and '#' not in vcf_line[0]:
                out_dict[vcf_line[0]].append(vcf_line)
    return out_dict

def split_by_chrom(vcf_filename):
    out_dict = split_by_chrom_helper(vcf_filename)
    # write each one to files
    for chrom, vcf_lines in out_dict.items():
        out_filename = f'{chrom}.pileup.tsv'
        with open(out_filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            csvwriter.writerows(vcf_lines)

from collections import defaultdict
import csv
import sys
argv = sys.argv
vcf_filename = argv[1]
split_by_chrom(vcf_filename)