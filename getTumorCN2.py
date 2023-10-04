#!/usr/bin/env python
# get copy number of tumor with a matched normal in the same vcf file
# use the read depth information in the matched normal to normalize read depth in the tumor sample
# Input VCF file format for this script:
col2idx = {'#CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3, 'QUAL': 4, 'normal_rd': 5, 'tumor_rd': 6, 'tumor_h1': 7, 'tumor_h2': 8, 'normal_h1': 9, 'normal_h2': 10}

class Variant(object):
    """This is an internal class object to store variant info."""
    def __init__(self, chrName, pos, nr_gt, nr_dep, tr_gt, tr_dep, th1_gt, th1_dep, th2_gt, th2_dep):
        self.chrName = chrName     # string
        self.pos     = int(pos)    # position of variant
        self.nr_gt   = int(nr_gt)  # 0:ref/ref, 1:ref/alt, 2:alt/alt
        self.nr_dep  = nr_dep      # list of length 2 with read depths of ref and alt allele in normal sample
        self.tr_gt   = int(tr_gt)  # 0:ref/ref, 1:ref/alt, 2:alt/alt
        self.tr_dep  = tr_dep      # list of length 2 with read depths of ref and alt allele in tumor sample
        self.th1_gt  = th1_gt
        self.th1_dep = th1_dep
        self.th2_gt  = th2_gt
        self.th2_dep = th2_dep
        
    def __str__(self):
        stout  = '\t'.join(['chrName', 'position','th1_gt', 'th1_dep', 'th2_gt', 'th2_dep'])
        stout += '\n'
        stout += '\t'.join([self.chrName,str(self.pos),str(self.th1_gt),str(self.th1_dep), str(self.th2_gt), str(self.th2_dep)])
        return stout

    def normalize_read_depth(self, normal_avg_depth):
        '''normalized tumor read depths using the read depth info from matched normal.
        process if the variant has ambiguous phasing'''
        if self.th1_gt == self.th2_gt and self.th1_dep == self.th2_dep:
            # print(self.chrName, self.pos)
            th1_dep_norm = -1
            th2_dep_norm = self.th2_dep/(sum(self.nr_dep)/normal_avg_depth)
        else:
            th1_dep_norm = self.th1_dep/(sum(self.nr_dep)/normal_avg_depth)
            th2_dep_norm = self.th2_dep/(sum(self.nr_dep)/normal_avg_depth)
        return th1_dep_norm, th2_dep_norm

def process_gt(gt):
    '''function to deal with gaps in genotypes:
    0: reference; 1: one ref/one alt; 2: alt; .: gap'''
    if '.' in gt: return '.'
    return sum([int(_) for _ in gt])

def gt_screening(read_gt, hap1_gt, hap2_gt, read_dep, depth_thres):
    r_gt  = process_gt(read_gt)
    h1_gt = process_gt(hap1_gt)
    h2_gt = process_gt(hap2_gt)
    if read_dep[0] <= depth_thres and r_gt == 1:
        r_gt = 2
    elif read_dep[1] <= depth_thres and r_gt == 1:
        r_gt = 0
    return r_gt, h1_gt, h2_gt

def get_info(parsed_line, depth_thres=3):
    '''Retrieve genotype and depth info from a row in vcf file format to a Variant object.
    For getHapDepth without a matched normal, the current parsed_line must contain 
    a heterozygous variant to get hap_depth.
    '''
    variant = None
    chrName = parsed_line[col2idx['#CHROM']]
    pos     = parsed_line[col2idx['POS']]
    normalr_gt  = parsed_line[col2idx['normal_rd']].split(':')[0].split('/')
    normalh1_gt = parsed_line[col2idx['normal_h1']].split(':')[0].split('/')
    normalh2_gt = parsed_line[col2idx['normal_h2']].split(':')[0].split('/')
    # if gap in normal sample, then it shouldn't be trusted
    if '.' in normalr_gt or '.' in normalh1_gt or '.' in normalh2_gt: return None
    # save info in matched normal sample
    nr_dep = [int(_) for _ in parsed_line[col2idx['normal_rd']].split(':')[1].split(',')]
    nr_gt, nh1_gt, nh2_gt = gt_screening(normalr_gt, normalh1_gt, normalh2_gt, nr_dep, depth_thres)
    # normal read-based genotype is different from normal assembly-based genotypes
    if nr_gt == 1 and nh1_gt == nh2_gt: return None
    if nh1_gt == 1 or nh2_gt == 1: return None
    
    # processing tumor samples
    tumorr_gt  = parsed_line[col2idx['tumor_rd']].split(':')[0].split('/')
    tumorh1_gt = parsed_line[col2idx['tumor_h1']].split(':')[0].split('/')
    tumorh2_gt = parsed_line[col2idx['tumor_h2']].split(':')[0].split('/')
    # if gap in tumor read or no coverage in either haplotype
    if '.' in tumorr_gt or ('.' in tumorh1_gt and '.' in tumorh2_gt): return None
    # save info in tumor sample
    tr_dep = [int(_) for _ in parsed_line[col2idx['tumor_rd']].split(':')[1].split(',')]
    tr_gt, th1_gt, th2_gt = gt_screening(tumorr_gt, tumorh1_gt, tumorh2_gt, tr_dep, depth_thres)
    # tumor read-based genotype is different from tumor assembly-based genotypes
    # **Note**: if one tumor chromosome has odd number of copy number, 
    # th1_gt or th2_gt might be 1 since one haplotype from h2 may go to h1...?
    if tr_gt == 1 and (th1_gt == th2_gt or th1_gt == '.' or th2_gt == '.'): return None
    if th1_gt == 1 or th2_gt == 1: return None

    # pass variant "QC", saving information...
    # "hemizygous" in tumor_hap2
    th1_dep = th2_dep = 0
    if th1_gt == '.':
        th1_dep = 0
        th2_dep = tr_dep[int(th2_gt/2)]
    # "hemizygous" in tumor_hap1
    elif th2_gt == '.':
        th2_dep = 0
        th1_dep = tr_dep[int(th1_gt/2)]
    # both tumor haps "homozygous": we can't determine the phase of these based on the hifiasm assembly
    # will be dealed later in function read_depth_from_vcf
    elif th1_gt == th2_gt:
        th1_dep = tr_dep[int(th1_gt/2)]
        th2_dep = tr_dep[int(th2_gt/2)]
    # good heterozygous variant in tumor sample!
    elif tr_gt == 1 and th1_gt != th2_gt:
        th1_dep = tr_dep[int(th1_gt/2)]
        th2_dep = tr_dep[int(th2_gt/2)]
    # else:
        # print(f"Something wrong with hap/get_info. Check {parsed_line}.")

    variant = Variant(chrName, pos, nr_gt, nr_dep, tr_gt, tr_dep, th1_gt, th1_dep, th2_gt, th2_dep)
    return variant

def read_depth_from_vcf(vcf_filename, normal_avg_depth, depth_thres=1):
    '''Returns a list of read depth, where each contains the position of variant, hap1 read depth, 
    and hap2 read depth. The read depths are normalized using the read depth info from matched normal.
    '''
    variants = []
    read_depth = []
    for line in open(vcf_filename, "r"):
        parsed_line = line.strip("\n").split("\t")
        if len(parsed_line) < 2:
            continue
        else:
            variants.append(get_info(parsed_line, depth_thres))

    for variant in variants:
        if variant:
            hap_depth = variant.normalize_read_depth(normal_avg_depth)
            if -1 in hap_depth: # ambiguous phasing!
                continue
            read_depth.append([variant.pos, *hap_depth])
    return read_depth

def smooth_allele_freq(read_depth, window_num, window_size, overlap_size=2):
    span = read_depth[-1][0]
    if window_num is not None:
        window = int(span/window_num)
    else:
        window = window_size
    res = []
    j = 0
    for i in range(window, span+window, int(window/overlap_size)):
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
            
    # median_depth = stats.median(total_depth)
    # print(min(total_depth))
    # print(f'median depth: {median_depth}')
    return res

def smooth_copy_number(read_depth, window_num, window_size, avg_depth, overlap_size):
    span = read_depth[-1][0]
    if window_num is not None:
        window = int(span/window_num)
    else:
        window = window_size
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
            # bug line: if j_flag and curr[0] > i+window/overlap_size:
            if j_flag and curr[0] >= i:
                next_j = j
                j_flag = False
                break
            elif curr[0] < i:
                hap1_depth += curr[1]
                hap2_depth += curr[2]
                j += 1
                occurence += 1
            else:
                print(f"Error in smooth copy number at position {curr[0]} in bin {i}.")
        if hap1_depth != 0 or hap2_depth != 0:
            depth_sum = hap1_depth + hap2_depth
            # avg_depth = depth_sum / occurence / 4
            res.append([i, hap1_depth/occurence/avg_depth, hap2_depth/occurence/avg_depth])
    return res

def get_CN_data(vcf_filename, tumor_avg_depth, normal_avg_depth, window_num=None, window_size=1000000, overlap_size=2):
    # fixed size 
    '''
    Without a matched normal, CN and AF are calculated from the position of heterozygous SNP, 
    hap1 read depth, and hap2 read depth from a pileup VCF file.
    Input:
        vcf_filename:   STR     a pileup VCF file of one chromosome from a sample (tumor/normal).
        avg_depth:      INT     the HiFi sequencing read depth divided by ploidy level.
        window_num:     INT     number of windows to bin the CN across one chromosome. Default is None.
                                If it is not none, then window_size is determined by window_num instead
        window_size:    INT     size of each bin/window across one chromosome. Default is 1Mb.
                                Note that this variable takes effect only when window_num is None. 
        overlap_size:   INT     size overlap for window sliding when calculating CN within a bin.
                                Default is 2, meaning the window slides 1/2 of bin size.
    Output:
        CN: a 3D list containing the index of each bin, the average CN of Hap1 and Hap2 in each bin, respectively.
    '''
    read_depth = read_depth_from_vcf(vcf_filename, normal_avg_depth)
    if window_num is not None:
        CN = smooth_copy_number(read_depth, window_num, None, tumor_avg_depth, overlap_size)
        AF = smooth_allele_freq(read_depth, window_num, None, overlap_size)
    else:
        CN = smooth_copy_number(read_depth, None, window_size, tumor_avg_depth, overlap_size)
        AF = smooth_allele_freq(read_depth, None, window_size, overlap_size)
    # CN = np.array(CN).T
    return list(CN), list(AF)

def get_tumorCN_matchedNormal(tumor_vcf_filename, tumor_avg_depth, normal_vcf_filename, normal_avg_depth, window_num=None, window_size=1000000, overlap_size=2):
    '''outdated
    Tumor sample CN calculated with a matched normal. 
    New Input:
        normal_read_depth: the output of read_depth_from_vcf from the matching normal VCF file
    '''
    tumor_read_depth, normal_read_depth = read_depth_matchedNormal(tumor_vcf_filename, normal_vcf_filename)
    if window_num is not None:
        tumor_CN  = smooth_copy_number(tumor_read_depth, window_num, None, tumor_avg_depth, overlap_size)
        tumor_AF  = smooth_allele_freq(tumor_read_depth, window_num, None, overlap_size)
        normal_CN = smooth_copy_number(normal_read_depth, window_num, None, normal_avg_depth, overlap_size)
        normal_AF = smooth_allele_freq(normal_read_depth, window_num, None, overlap_size)
    else:
        tumor_CN  = smooth_copy_number(tumor_read_depth, None, window_size, tumor_avg_depth, overlap_size)
        tumor_AF  = smooth_allele_freq(tumor_read_depth, None, window_size, overlap_size)
        normal_CN = smooth_copy_number(normal_read_depth, None, window_size, normal_avg_depth, overlap_size)
        normal_AF = smooth_allele_freq(normal_read_depth, None, window_size, overlap_size)
    # tumor_CN = np.array(tumor_CN).T
    # normal_CN = np.array(normal_CN).T
    output = [tumor_CN, tumor_AF, normal_CN, normal_AF]
    return output

def write2csv(alist, out_filename):
    with open(out_filename, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(alist)

def get_chrName(filename):
    tmp = filename.split('/')[-1].split('.')
    chrom = ['chr' in i for i in tmp]
    chrName = [tmp[i] for i in range(len(chrom)) if chrom[i]][0]
    return chrName

def usage():
    print("Usage: getTumorCN.py --tumor_read_depth [HiFi read depth of tumor sample] -t [directory to tumor pileup_byChr files]")
    print("Options:")
    print("  -w INT      Number of windows for each chromosome")
    print("  -s INT      Number of sliding windows for each bin")
    print("  -t STR      Directory containing pileup VCF files for one chromosome of tumor sample")
    print("  -n STR      Directory containing pileup VCF files for one chromosome of matched normal sample")
    print("  --tumor-read-depth     INT      HiFi read depth of tumor sample")
    print("  --normal-read-depth    INT      HiFi read depth of matched normal sample")
    print("Note: processed pileup VCF file split by chromosomes should have the format <chrName>.pileup.tsv")
    print("Provide matched normal sample to resolve CNV in hemizygous regions.")

def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:],"f:w:s:", 
                               ["tumor-read-depth=","normal-read-depth="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) < 1:
        usage()
        sys.exit(1)
    
    vcf_filenames  = []
    window_num   = None
    overlap_size = 2
    for opt, arg in opts:
        if opt == '-f':
            for filename in glob(f'{arg}/chr*.pileup.tsv'):
                vcf_filenames.append(filename)
        elif opt == '-w': 
            print('Dynamic window size instead of static (1Mb)...')
            print(f'Number of windows for each chromosome: {arg}')
            window_num   = int(arg)
        elif opt == '-s': 
            print(f'Number of sliding windows for each bin: {arg}')
            overlap_size = int(arg)
        elif opt in ['--tumor-read-depth']:  tumor_avg_depth = float(arg)
        elif opt in ['--normal-read-depth']: normal_avg_depth = float(arg)
    
    if len(vcf_filenames) == 0:
        print("VCF file error.")
        sys.exit(1)

    else:
        for filename in vcf_filenames:
            print(filename)
            chrName = get_chrName(filename)
            # getCN&AF
            tumor_CN,tumor_AF = get_CN_data(filename, tumor_avg_depth, normal_avg_depth, window_num=window_num, overlap_size=overlap_size)
            write2csv(tumor_CN, f'{chrName}.tumor.CN.tsv')
            write2csv(tumor_AF, f'{chrName}.tumor.AF.tsv')


if __name__ == "__main__":
    import os
    import sys
    import csv
    import getopt
    from glob import glob
    # from statistics import mean
    # import matplotlib.pyplot as plt
    from collections import defaultdict
    main(sys.argv)