# take out the copy number part from the older version and only include heterozygosity evaluation
import numpy as np
import os
import sys
import csv
from statistics import mean
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib import ticker
from matplotlib.patches import Rectangle
import matplotlib as mpl
params = {'axes.spines.top': False, 'axes.spines.right': False}
mpl.rcParams.update(params)

# from preprocess_vcf.ipynb
col2idx = {'#CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3, 'QUAL': 4, 'read': 5, 'hap1': 6, 'hap2': 7}


# ------ Section 1: get allelic/genotypic information of each variant across the genome ------
class Variant(object):
    """This is an internal class object to store variant info."""
    def __init__(self, chrName, pos, read, hap1, hap2):
        self.chrName  = chrName  # string
        self.pos      = int(pos)
        self.read     = int(read)     
        self.hap1     = hap1
        self.hap2     = hap2
        
    def __str__(self):
        return '\t'.join([self.chrName,str(self.pos),str(self.read),str(self.hap1),str(self.hap2)])

def get_gt(sample_info, get_depth = False):
    '''Retrieve genotype info from the sample column in vcf file format to a list'''
    tmp = sample_info.split(':')[0].split('/')
    if '.' in tmp:
        gt = '.'
    else:
        gt = sum([int(_) for _ in tmp])
    if get_depth:
        return gt, [int(_) for _ in sample_info.split(':')[1].split(',')]
    return gt

def read_genotype_from_vcf(file_name, hm_thres = 3, out_file = None):
    '''Divide homozygous/heterozygous reads from processed vcf file. 
       If one of the heterozygous genotype has read depth smaller than the threshold, 
       it is considered as a homozygous position.
       This function can also combine both homotypes and heterotypes and write to a tsv file.
       
    Input:
        - file_name: vcf file input name
        - hm_thres: Filtering threshold of read depth for ambiguous homotype. Default is 3.
        
    Output:
        - homotypes: list of variants where the reads are homozygous -- carry 0 or 2 ALT alleles
        - heterotypes: list of variants where the reads are heterozygous -- carry 1 ALT alleles
        - out_file: genotypes output file with homotypes and heterotypes combined
    '''
    print(f"Running read_genotype({file_name})...")
    homotypes = []
    heterotypes = []
    for line in open(file_name, "r"):
        parsed_line = line.strip("\n").split("\t")
        if len(parsed_line) < 2:
            continue
        else:
            rd_gt, depths = get_gt(parsed_line[col2idx['read']], True)
            hap1_gt = get_gt(parsed_line[col2idx['hap1']])
            hap2_gt = get_gt(parsed_line[col2idx['hap2']])
            
            # 1) homozygous if rd_gt is 0/0 or 1/1
            if rd_gt == 0 or rd_gt == 2:
                variant = Variant(*parsed_line[0:2], rd_gt, hap1_gt, hap2_gt)
                homotypes.append(variant)
            # 2) homozygous if one of read depths is smaller than threshold
            elif depths[0] < 3:
                rd_gt = 2
                variant = Variant(*parsed_line[0:2], rd_gt, hap1_gt, hap2_gt)
                homotypes.append(variant)
            elif depths[1] < 3:
                rd_gt = 0
                variant = Variant(*parsed_line[0:2], rd_gt, hap1_gt, hap2_gt)
                homotypes.append(variant)
            # 3) heterozygous otherwise
            else:
                variant = Variant(*parsed_line[0:2], rd_gt, hap1_gt, hap2_gt)
                heterotypes.append(variant)
            if out_file:
                out_file.write(variant.__str__() + '\n')
    print(f"Number of read homotypes: {len(homotypes)}.")
    print(f"Number of read heterotypes: {len(heterotypes)}.")
    return homotypes, heterotypes

def hm_analysis(homotypes):
    '''
    Cases:
        1. Both haps captured the homozygosity
        2. One hap captured, one is '.'
        2. Clustering gaps
        3. False duplication
        4. Missing homozygosity
    '''
    # stored as dictionary of lists of positions with chromosome as keys
    hm_dicap  = defaultdict(list)
    hm_unicap = defaultdict(list)
    hm_miss   = defaultdict(list)
    hm_gap    = defaultdict(list)
    hm_false  = defaultdict(list)
    count = defaultdict(int)
    
    for variant in homotypes:
        count[variant.chrName] += 1
        if variant.read == variant.hap1 == variant.hap2:
            hm_dicap[variant.chrName].append(variant.pos)
        elif variant.hap1 == '.' and variant.hap2 == '.':
            hm_gap[variant.chrName].append(variant.pos)
        elif variant.hap1 == 1 or variant.hap2 == 1:
            hm_false[variant.chrName].append(variant.pos)
        elif variant.hap1 == '.' or variant.hap2 == '.':
            hm_unicap[variant.chrName].append(variant.pos)
        else:
            hm_miss[variant.chrName].append(variant.pos)
    return hm_dicap, hm_unicap, hm_miss, hm_gap, hm_false, count

def ht_analysis(heterotypes):
    '''
    Cases:
        1. Both haps captured the heterozygosity
        2. Clustering gaps
        3. False duplication
        4. Missing heterozygosity
    '''
    # stored as dictionary of lists of positions with chromosome as keys
    ht_true  = defaultdict(list) 
    ht_gap   = defaultdict(list)
    ht_false = defaultdict(list)
    ht_miss = defaultdict(list)
    count = defaultdict(int)
    
    for variant in heterotypes:
        count[variant.chrName] += 1
        if (variant.hap1 == 2 and variant.hap2 == 0) or (variant.hap1 == 0 and variant.hap2 == 2):
            ht_true[variant.chrName].append(variant.pos)
        elif variant.hap1 == '.' or variant.hap2 == '.':
            ht_gap[variant.chrName].append(variant.pos)
        elif variant.hap1 == 1 or variant.hap2 == 1:
            ht_false[variant.chrName].append(variant.pos)
        else:
            ht_miss[variant.chrName].append(variant.pos)
    return ht_true, ht_gap, ht_false, ht_miss, count

def read_ctgaln_from_bed(bed_filename):
    '''Bed file stores info about all chromosome coverages from the assembled contigs
    '''
    chr_cvg = defaultdict(list)
    for line in open(bed_filename,'r'):
        parsed_line = line.strip("\n").split("\t")
        if len(parsed_line) < 2:
            continue
        else:
            chr_cvg[parsed_line[0]].append([int(parsed_line[1]), int(parsed_line[2])])
    return chr_cvg

def get_frac_within_window(curr_cases, start, end):
    '''
    A helper method for cnt2frac function.
    '''
    cnt  = defaultdict(float)
    frac = []
    for i in range(len(curr_cases)):
        case = curr_cases[i]
        for pos in case:
            if pos >= start and pos < end:
                cnt[i] += 1
        cnt[i] += sys.float_info.epsilon
    total_cnt = sum(cnt.values())
    if total_cnt < 1:
        return [0 for i in range(len(curr_cases))]
    for c in cnt.values():
        frac.append(c/total_cnt)
    return frac

def cnt2frac(hmt_all, window = 1000000):
    '''Older version which treat all chromosome with different length all with the same window size.
    The input hmt_all is a list of cases with ascending-ordered positions of either homotypes or heterotypes.
    We want to count all homozygote/heterozygote cases within that window (default 1mb), 
    and have the fraction of each case.
    '''
    chroms = [i for i in hmt_all[0].keys() if '_' not in i]
    output = defaultdict(list)
    
    for chrom in chroms:
        curr_cases = [hmt_case[chrom] for hmt_case in hmt_all]
        curr_maxlen = max([max(case) for case in curr_cases if len(case)>0])
        start = 0
        end = start + window + 1
        
        while start < curr_maxlen:
            frac = get_frac_within_window(curr_cases, start, end)
            output[chrom].append(frac)
            start += window
            end = start + window + 1
    return output

def cnt2frac_byChr(hmt_all, chrom, window_num = 250):
    '''
    The input hmt_all is a list of cases with ascending-ordered positions of either homotypes or heterotypes.
    We want to count all homozygote/heterozygote cases within a window size 
    such that we have a default of 250 windows, and have the fraction of each case.
    '''
    output = []
    
    curr_cases = [hmt_case[chrom] for hmt_case in hmt_all]
    curr_maxlen = max([max(case) for case in curr_cases if len(case)>0])
    window = int(curr_maxlen/window_num)
    start = 0
    end = start + window + 1

    while start < curr_maxlen:
        frac = get_frac_within_window(curr_cases, start, end)
        output.append(frac)
        start += window
        end = start + window + 1
    
    # use this returned window size for plotting
    return np.array(output).T, window



# ------ Section 4: Plotting ------

def integrate_plotting_helper(hmt_fraction, window, chrom):
    # Fig size
    fig, ax = plt.subplots(figsize=(15, 4))
    # hmt_plotting
    xt = np.arange(0, len(hmt_fraction[0])*window, window)
    divide_array = np.array([1000000]*len(xt))
    x = np.divide(xt,divide_array).tolist()
    window = window/1000000
#     print(hmt_fraction[0])
    ax.bar(x, hmt_fraction[0], width=window, align='edge', color='blue', alpha=0.7) # hm_dicap
    ax.bar(x, hmt_fraction[1], width=window, bottom=hmt_fraction[0], align='edge', color='green', alpha=0.5) # hm_unicap
    ax.bar(x, hmt_fraction[2], width=window, bottom=hmt_fraction[0]+hmt_fraction[1], align='edge', color='purple', alpha=0.5) # hm_false
    ax.bar(x, hmt_fraction[3], width=window, bottom=hmt_fraction[0]+hmt_fraction[1]+hmt_fraction[2], align='edge', color='red', alpha=0.5) # hm_miss
    ax.bar(x, hmt_fraction[4], width=window, bottom=hmt_fraction[0]+hmt_fraction[1]+hmt_fraction[2]+hmt_fraction[3], align='edge', color='silver', alpha=0.5) # hm_gap
    ax.bar(x, hmt_fraction[5]*-1, width=window, align='edge', color='blue',alpha=0.7) # ht_true
    ax.bar(x, hmt_fraction[6]*-1, width=window, bottom=hmt_fraction[5]*-1, align='edge', color='purple', alpha=0.5) # ht_false
    ax.bar(x, hmt_fraction[7]*-1, width=window, bottom=(hmt_fraction[5]+hmt_fraction[6])*-1, align='edge', color='red', alpha=0.5) # ht_miss
    ax.bar(x, hmt_fraction[8]*-1, width=window, bottom=(hmt_fraction[5]+hmt_fraction[6]+hmt_fraction[7])*-1, align='edge', color='silver', alpha=0.5) # ht_gap
    
    plt.ylabel("Ratio")
    plt.legend([Rectangle((0,0),3,1.5, facecolor='blue',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='green',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='purple',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='red',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='silver',edgecolor='black')],
                     ['di', 'uni', 'het', 'non', 'gap'])

    xlim = x[-1]+10
    plt.xlabel('Position (Mb)')
    plt.xticks(np.arange(0,xlim,10))
    plt.ticklabel_format(style='plain')
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    plt.savefig(f'hetEval.{chrom}.pdf',format='pdf',bbox_inches='tight')
    plt.close()

def integrate_plotting(hmt_all, chrom, window_num=250):
    hmt_fraction, window = cnt2frac_byChr(hmt_all, chrom, window_num)
    integrate_plotting_helper(hmt_fraction, window, chrom)
    return None


# Everything together
def haplotype_eval(processed_vcf):
    # check if pileup_byChr exists. If not, create directory and split vcf files.
    homotypes, heterotypes = read_genotype_from_vcf(processed_vcf)
    hm_dicap, hm_unicap, hm_miss, hm_gap, hm_false, hm_count = hm_analysis(homotypes)
    ht_true, ht_gap, ht_false, ht_miss, ht_count = ht_analysis(heterotypes)
    hmt_all = [hm_dicap, hm_unicap, hm_false, hm_miss, hm_gap, 
               ht_true, ht_false, ht_miss, ht_gap]
    hmt_idx2col = {'hm_dicap':0, 'hm_unicap':1, 'hm_false':2, 'hm_miss':3, 'hm_gap':4, 
                   'ht_true':5, 'ht_false':6, 'ht_miss':7, 'ht_gap':8}
    chroms = [i for i in hmt_all[0].keys() if '_' not in i and 'M' not in i]
    
    for chrom in chroms:
        print(chrom)
        integrate_plotting(hmt_all, chrom, window_num=200)
        
        
argv = sys.argv
if len(argv) < 2 or len(argv) > 3:
    print('Usage: haplotypeEval.py [processed_vcf]')
    sys.exit(1)
print(f'Running {argv[0]} on file {argv[1]}')
haplotype_eval(argv[1])

