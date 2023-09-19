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

# ------ Section 2: get copy number information from each variant across the genome ------
def split_by_chrom_helper(vcf_filename):
    '''helper to split vcf file rows by chromosomes'''
    out_dict = defaultdict(list)
    with open(vcf_filename,'r') as vcf_file:
        for line in vcf_file.readlines():
            vcf_line = line.rstrip().split('\t')
            if '_' not in vcf_line[0]:
                out_dict[vcf_line[0]].append(vcf_line)
    return out_dict

def split_by_chrom(vcf_filename, directory = '.'):
    out_dict = split_by_chrom_helper(vcf_filename)
    # write each one to files
    for chrom, vcf_lines in out_dict.items():
        out_filename = f'{directory}/{chrom}.pileup.tsv'
        with open(out_filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            csvwriter.writerows(vcf_lines)

def get_each_hap_depth(parsed_line, depth_thres=3):
    '''Retrieve genotype and depth info from a row in vcf file format to a list'''
    read_gt = parsed_line[col2idx['read']].split(':')[0].split('/')
    hap1_gt = parsed_line[col2idx['hap1']].split(':')[0].split('/')
    hap2_gt = parsed_line[col2idx['hap2']].split(':')[0].split('/')
    
    # either gap or homozygotes or overlapped contig
    if '.' in read_gt or '.' in hap1_gt or '.' in hap2_gt:
        return None
    
    rd_gt =  sum([int(_) for _ in read_gt])
    h1_gt =  sum([int(_) for _ in hap1_gt])
    h2_gt =  sum([int(_) for _ in hap2_gt])
    read_depth = [int(_) for _ in parsed_line[col2idx['read']].split(':')[1].split(',')]
    
    if rd_gt == 0 or rd_gt == 2:
        return None
    elif read_depth[0] <= depth_thres or read_depth[1] <= depth_thres:
        return None
    # **Note**: 
    # if one cancer chromosome has odd number of copy number, 
    # h1_gt or h2_gt might be 1 since one haplotype from h2 may go to h1...?
    elif h1_gt == 1 or h2_gt == 1:
        return None
    elif h1_gt == h2_gt:
        # print(f"Missing heterozygosity at {parsed_line[col2idx['#CHROM']]}: {parsed_line[col2idx['POS']]}")
        return None
    else:
        h1_depth = read_depth[int(h1_gt/2)]
        h2_depth = read_depth[int(h2_gt/2)]
        return [h1_depth, h2_depth]

def read_depth_from_vcf(file_name):
    '''Returns a list of read depth. 
    Each element contains the position of heterozygous SNP, hap1 read depth, and hap2 read depth.'''
    read_depth = []
    for line in open(file_name, "r"):
        parsed_line = line.strip("\n").split("\t")
        if len(parsed_line) < 2:
            continue
        else:
            hap_depth = get_each_hap_depth(parsed_line)
            if hap_depth:
                read_depth.append([int(parsed_line[col2idx['POS']]), *hap_depth])
    return read_depth

def smooth_allele_freq(read_depth, window=1000000):
    span = read_depth[-1][0]
    res  = []
    j = 0
    for i in range(window, span+window, window):
        hap1_depth = hap2_depth = 0
        while j < len(read_depth):
            curr = read_depth[j]
            if curr[0] < i:
                hap1_depth += curr[1]
                hap2_depth += curr[2]
                j += 1
            else:
                break
        if hap1_depth != 0:
            depth_sum = hap1_depth + hap2_depth
            res.append([i, hap1_depth/depth_sum, hap2_depth/depth_sum])
            
    median_depth = stats.median(total_depth)
#     print(min(total_depth))
#     print(f'median depth: {median_depth}')
    return res

def plot_allele_freq(file_name, window=1000000):
    tmp = file_name.split('/')[-1].split('.')
    chrom = ['chr' in i for i in tmp]
    chrName = [tmp[i] for i in range(len(chrom)) if chrom[i]][0]

    read_depth = read_depth_from_vcf(file_name)
    freq,median_depth = smooth_allele_freq(read_depth, window)
    freq = np.array(freq).T
    plt.figure(figsize=(20,5))
    plt.plot(freq[0], freq[1], linewidth=2,linestyle='dashed',label ='Hap1')
    plt.plot(freq[0], freq[2], linewidth=2,linestyle='dashed',label ='Hap2')

    plt.xlabel("Position")
    plt.ylabel("Allele Ratio")
    plt.xticks(np.arange(0,freq[0][-1]+1000000,10000000),rotation=45)
    plt.ticklabel_format(style='plain')
    plt.legend()
    plt.savefig(f'{chrName}.alleleFreq.pdf', format='pdf')
    plt.show()
    
    return list(freq)

def smooth_copy_number(read_depth, window_num, avg_depth, overlap_size):
    span = read_depth[-1][0]
    window = int(span/window_num)
    res  = []
    j = occurence = 0
    # In order to create a sliding window
    next_j = 0
    j_flag = True
    
    # window/overlap_size for overlapping sliding
    for i in range(window, span+window, int(window/overlap_size)):
        occurence = 0
        hap1_depth = hap2_depth = 0
        j = next_j
        j_flag = True
        while j < len(read_depth):
            curr = read_depth[j]
            # Record the next starting j
            if j_flag and curr[0] > i+window/overlap_size:
                next_j = j
                j_flag = False
                
            if curr[0] < i:
                hap1_depth += curr[1]
                hap2_depth += curr[2]
                j += 1
                occurence += 1
            else:
                break
        if hap1_depth != 0:
            depth_sum = hap1_depth + hap2_depth
            res.append([i, hap1_depth/occurence/avg_depth, hap2_depth/occurence/avg_depth])
            
    return res

def get_CN_data(file_name, avg_depth = 69/4, window_num = 250, overlap_size=2):
    # if only chromosome name is provided
    # the file must be in a subdirectory of the current script in a format of {chrom}.pileup.tsv
    directory = "pileup_byChr"
    if '.tsv' not in file_name and 'chr' in file_name:
        file_name = f'{directory}/{file_name}.pileup.tsv'

    tmp = file_name.split('/')[-1].split('.')
    chrom = ['chr' in i for i in tmp]
    chrName = [tmp[i] for i in range(len(chrom)) if chrom[i]][0]

    read_depth = read_depth_from_vcf(file_name)
    CN = smooth_copy_number(read_depth, window_num, avg_depth, overlap_size)
    CN = np.array(CN).T
    
    return list(CN)

# ------ Section 3: get contig length information from processed BED file from PAF across the genome ------
def get_ctginfo(cvg_list, window=10000):
    '''Get average contig length of one chromosome across its positions.
    Note that this length might be a partial contig block that's aligned to a part of the chromosome.'''
    ctg_len_tmp = defaultdict(list)
    ctg_len = dict()
    for cvg_pos in cvg_list:
        curr_ctg_len = cvg_pos[1] - cvg_pos[0]
        appr_start = cvg_pos[0] // window * window
        appr_end   = cvg_pos[1] // window * window
        for i in range(appr_start, appr_end, window):
            ctg_len_tmp[i].append(curr_ctg_len)
    for key, value in ctg_len_tmp.items():
        ctg_len[key] = mean(ctg_len_tmp[key])
    return ctg_len



# ------ Section 4: Plotting ------
def get_plotting_data(hmt_all, chr_cvg1, chr_cvg2, chrom, hifi_read_depth, window_num=250, ctg_window = 10000):
    '''get ctg_len, CN, and hmt_fraction for one designated chromosome'''
    hmt_fraction, window = cnt2frac_byChr(hmt_all, chrom, window_num)
    CN = get_CN_data(chrom, avg_depth=hifi_read_depth, window_num=window_num)
    cvg1 = chr_cvg1[chrom]
    cvg2 = chr_cvg2[chrom]
    ctg_len = (get_ctginfo(cvg1, window=ctg_window),get_ctginfo(cvg2, window=ctg_window))
    return ctg_len, CN, hmt_fraction, window, chrom

def integrate_plotting_helper(ctg_len, CN, hmt_fraction, window, chrom):
    # Fig size
    plt.rcParams["figure.figsize"]=16,12
    fig, axes = plt.subplots(nrows=3, sharex=True)
    
    # contig length
    ctg_len1, ctg_len2 = ctg_len
    axes[0].bar(ctg_len1.keys(), ctg_len1.values(), width=window, color='lavender')
    ctg_len2_key, ctg_len2_value = list(ctg_len2.keys()), list(ctg_len2.values())
    ctg_len2_value = np.array(ctg_len2_value) * -1
    axes[0].bar(ctg_len2_key, ctg_len2_value, width=window, color='thistle')
    
    axes[0].legend([Rectangle((0,0),3,1.5, facecolor='lavender',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='thistle',edgecolor='black')],
                     ['hap1', 'hap2'])
    axes[0].set_ylabel('Average contig alignment length')

    # Copy number
    axes[1].scatter(CN[0], CN[1], c='goldenrod', s=10, alpha=.7, label ='Hap1')
    axes[1].scatter(CN[0], CN[2], c='#5f7ba8', s=10, alpha=.7, label ='Hap2')
    axes[1].set_ylabel("Copy Number")
    axes[1].legend()
    
    # hmt_plotting
    x = np.arange(0, len(hmt_fraction[0])*window, window)
    
    axes[2].bar(x, hmt_fraction[0], width=window, align='edge', color='blue', alpha=0.7) # hm_dicap
    axes[2].bar(x, hmt_fraction[1], width=window, bottom=hmt_fraction[0], align='edge', color='green', alpha=0.5) # hm_unicap
    axes[2].bar(x, hmt_fraction[2], width=window, bottom=hmt_fraction[0]+hmt_fraction[1], align='edge', color='purple', alpha=0.5) # hm_false
    axes[2].bar(x, hmt_fraction[3], width=window, bottom=hmt_fraction[0]+hmt_fraction[1]+hmt_fraction[2], align='edge', color='red', alpha=0.5) # hm_miss
    axes[2].bar(x, hmt_fraction[4], width=window, bottom=hmt_fraction[0]+hmt_fraction[1]+hmt_fraction[2]+hmt_fraction[3], align='edge', color='silver', alpha=0.5) # hm_gap
    axes[2].bar(x, hmt_fraction[5]*-1, width=window, align='edge', color='blue',alpha=0.7) # ht_true
    axes[2].bar(x, hmt_fraction[6]*-1, width=window, bottom=hmt_fraction[5]*-1, align='edge', color='purple', alpha=0.5) # ht_false
    axes[2].bar(x, hmt_fraction[7]*-1, width=window, bottom=(hmt_fraction[5]+hmt_fraction[6])*-1, align='edge', color='red', alpha=0.5) # ht_miss
    axes[2].bar(x, hmt_fraction[8]*-1, width=window, bottom=(hmt_fraction[5]+hmt_fraction[6]+hmt_fraction[7])*-1, align='edge', color='silver', alpha=0.5) # ht_gap
    
    axes[2].set_ylabel("Ratio")
    axes[2].legend([Rectangle((0,0),3,1.5, facecolor='blue',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='green',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='purple',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='red',edgecolor='black'),
                     Rectangle((0,0),3,1.5, facecolor='silver',edgecolor='black')],
                     ['dicover', 'unicover', 'het_hap', 'noncover', 'gap_hap'])

    xlim = x[-1]+window+1
    plt.xlabel('Chromosome Coordinates')
    plt.xticks(np.arange(0,xlim,10000000),rotation=45)
    plt.ticklabel_format(style='plain')
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    plt.savefig(f'hmt.{chrom}.pdf',format='pdf',bbox_inches='tight')
    plt.close()

def integrate_plotting(hmt_all, chr_cvg1, chr_cvg2, chrom, hifi_read_depth, window_num=250, ctg_window = 10000):
    ctg_len, CN, hmt_fraction, window, chrom = get_plotting_data(hmt_all, chr_cvg1, chr_cvg2, 
                                                                 chrom, window_num, ctg_window)
    integrate_plotting_helper(ctg_len, CN, hmt_fraction, window, chrom)
    return None


# Everything together
def haplotype_eval(hifi_read_depth, processed_vcf, hap1aln_bed, hap2aln_bed):
    # check if pileup_byChr exists. If not, create directory and split vcf files.
    parent_dir = os.getcwd()
    directory = "pileup_byChr"
    path = os.path.join(parent_dir, directory)
    if not os.path.exists(path):
        os.mkdir(path)
        split_by_chrom(processed_vcf, directory=directory)
    elif len(os.listdir(path)) == 0:
        split_by_chrom(processed_vcf, directory=directory)
        
    homotypes, heterotypes = read_genotype_from_vcf(processed_vcf)
    chr_cvg1 = read_ctgaln_from_bed(hap1aln_bed)
    chr_cvg2 = read_ctgaln_from_bed(hap2aln_bed)
    hm_dicap, hm_unicap, hm_miss, hm_gap, hm_false, hm_count = hm_analysis(homotypes)
    ht_true, ht_gap, ht_false, ht_miss, ht_count = ht_analysis(heterotypes)
    hmt_all = [hm_dicap, hm_unicap, hm_false, hm_miss, hm_gap, 
               ht_true, ht_false, ht_miss, ht_gap]
    hmt_idx2col = {'hm_dicap':0, 'hm_unicap':1, 'hm_false':2, 'hm_miss':3, 'hm_gap':4, 
                   'ht_true':5, 'ht_false':6, 'ht_miss':7, 'ht_gap':8}
    chroms = [i for i in hmt_all[0].keys() if '_' not in i and 'M' not in i]
    
    for chrom in chroms:
        print(chrom)
        integrate_plotting(hmt_all, chr_cvg1, chr_cvg2, chrom, hifi_read_depth, window_num=200, ctg_window = 10000)
        
        
argv = sys.argv
if len(argv) < 5 or len(argv) > 6:
    print('Usage: haplotypeEval.py [HiFi_read_depth] [processed_vcf] [hap1aln_bed] [hap2aln_bed]')
    sys.exit(1)
print(f'Running {argv[0]} on file(s) {argv[2:5]}')
haplotype_eval(*argv[1:5])

