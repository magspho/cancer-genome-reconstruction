#!/usr/bin/env python
# get copy number of tumor from matched normal
col2idx = {'#CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3, 'QUAL': 4, 'read': 5, 'hap1': 6, 'hap2': 7}

def get_hap_depth(parsed_line, depth_thres=3):
    '''Retrieve genotype and depth info from a row in vcf file format to a list.
    For getHapDepth without a matched normal, the current parsed_line must contain 
    a heterozygous variant to get hap_depth. 
    '''
    read_gt = parsed_line[col2idx['read']].split(':')[0].split('/')
    hap1_gt = parsed_line[col2idx['hap1']].split(':')[0].split('/')
    hap2_gt = parsed_line[col2idx['hap2']].split(':')[0].split('/')
    
    # either gap or homozygotes or overlapped contig
    if '.' in read_gt or '.' in hap1_gt or '.' in hap2_gt:
        return None
    rd_gt =  process_gt(read_gt)
    h1_gt =  process_gt(hap1_gt)
    h2_gt =  process_gt(hap2_gt)
    read_depth = [int(_) for _ in parsed_line[col2idx['read']].split(':')[1].split(',')]
    
    if rd_gt == 2 or read_depth[0] <= depth_thres or read_depth[1] <= depth_thres:
        return None
    # **Note**: if one cancer chromosome has odd number of copy number, 
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

def process_gt(gt):
        '''function to deal with gaps in genotypes:
        0: reference; 1: one ref/one alt; 2: alt; .: gap
        '''
        if '.' in gt: return '.'
        return sum([int(_) for _ in gt])

def get_hap_depth_matchedNormal(parsed_line, depth_thres=3):
    '''If we have a matched normal (mode_tumor=True), then we find the correct 
    hifiasm-assembled tumor haplotype with the help of corresponding heterozygous variant 
    from the matched normal even in case of hemizygosity in tumor sample.'''
    read_gt = parsed_line[col2idx['read']].split(':')[0].split('/')
    hap1_gt = parsed_line[col2idx['hap1']].split(':')[0].split('/')
    hap2_gt = parsed_line[col2idx['hap2']].split(':')[0].split('/')
    # gap in read
    if '.' in read_gt or ('.' in hap1_gt and '.' in hap2_gt):
        return None
    rd_gt =  process_gt(read_gt)
    h1_gt =  process_gt(hap1_gt)
    h2_gt =  process_gt(hap2_gt)
    read_depth = [int(_) for _ in parsed_line[col2idx['read']].split(':')[1].split(',')]

    if h1_gt == '.' or h2_gt == '.' or (h1_gt == 2 and h2_gt == 2):
        # hemizygous
        if rd_gt == 2 or read_depth[0] <= depth_thres or read_depth[1] <= depth_thres:
            return [h1_gt, h2_gt, sum(read_depth)]
        else: return None
    else:
    # **Note**: if one cancer chromosome has odd number of copy number, 
    # h1_gt or h2_gt might be 1 since one haplotype from h2 may go to h1...?
        if h1_gt == 1 or h2_gt == 1:
            return None
        elif h1_gt == h2_gt:
            # print(f"Missing heterozygosity at {parsed_line[col2idx['#CHROM']]}: {parsed_line[col2idx['POS']]}")
            return None
        else:
            h1_depth = read_depth[int(h1_gt/2)]
            h2_depth = read_depth[int(h2_gt/2)]
            return [h1_depth, h2_depth]

def read_depth_matchedNormal(tumor_vcf_filename, normal_vcf_filename, depth_thres=3):
    '''Note the variant could be heterozygous or hemizygous in tumor sample with a matched normal.
    '''
    normal_heterotypes = dict()
    normal_read_depth  = []
    
    with open(normal_vcf_filename,'r') as normal_vcf:
        for line in normal_vcf.readlines():
            parsed_line = line.strip("\n").split("\t")
            if len(parsed_line) < 2:
                continue
            hap_depth = get_hap_depth(parsed_line)
            if hap_depth:
                curr_pos = int(parsed_line[col2idx['POS']])
                normal_read_depth.append([curr_pos, *hap_depth])
                normal_heterotypes[curr_pos] = parsed_line
    tumor_read_depth = []    
    with open(tumor_vcf_filename,'r') as tumor_vcf:
        for line in tumor_vcf.readlines():
            parsed_line = line.strip("\n").split("\t")
            if len(parsed_line) < 2:
                continue
            curr_pos = int(parsed_line[col2idx['POS']])
            hap_depth = get_hap_depth_matchedNormal(parsed_line)
            if hap_depth is None:
                continue
            # hemizygous and only one hap covers
            if '.' in hap_depth:
                if curr_pos in normal_heterotypes:
                    nongap_hap = 1 - hap_depth.index('.') # find which hap doesn't have gap
                    if hap_depth[nongap_hap] != 2:
                        print(f"Something wrong with hap. Check {parsed_line}.")
                        continue
                    # point read depth to that nongap_hap
                    else:
                        if nongap_hap == 0: # if hap1 doesn't have a gap
                            tumor_read_depth.append([curr_pos, hap_depth[2], 0])
                        else:
                            tumor_read_depth.append([curr_pos, 0, hap_depth[2]])
                else:
                    # the current tumor hemizygous variant is not covered in the matched normal
                    continue
            elif hap_depth[0] == hap_depth[1] == 2:
                tumor_read_depth.append([curr_pos, hap_depth[2]/2, hap_depth[2]/2])
            elif len(hap_depth) == 2:
                tumor_read_depth.append([curr_pos, *hap_depth])

    return tumor_read_depth, normal_read_depth

def read_depth_from_vcf(vcf_filename):
    '''Returns a list of read depth. 
    Each element contains the position of heterozygous SNP, hap1 read depth, and hap2 read depth.'''
    read_depth = []
    for line in open(vcf_filename, "r"):
        parsed_line = line.strip("\n").split("\t")
        if len(parsed_line) < 2:
            continue
        else:
            hap_depth = get_hap_depth(parsed_line)
            if hap_depth:
                read_depth.append([int(parsed_line[col2idx['POS']]), *hap_depth])
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
#     print(min(total_depth))
#     print(f'median depth: {median_depth}')
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

def get_CN_data(vcf_filename, avg_depth, window_num=None, window_size=1000000, overlap_size=2):
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
    read_depth = read_depth_from_vcf(vcf_filename)
    if window_num is not None:
        CN = smooth_copy_number(read_depth, window_num, None, avg_depth, overlap_size)
        AF = smooth_allele_freq(read_depth, window_num, None, overlap_size)
    else:
        CN = smooth_copy_number(read_depth, None, window_size, avg_depth, overlap_size)
        AF = smooth_allele_freq(read_depth, None, window_size, overlap_size)
    # CN = np.array(CN).T
    return list(CN), list(AF)

def get_tumorCN_matchedNormal(tumor_vcf_filename, tumor_avg_depth, normal_vcf_filename, normal_avg_depth, window_num=None, window_size=1000000, overlap_size=2):
    '''
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

def writeCN2file(CN, out_filename):
    with open(out_filename, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(CN)

def writeAF2file(AF, out_filename):
    with open(out_filename, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(AF)

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
        opts, args = getopt.getopt(argv[1:],"t:n:w:s:", 
                               ["tumor-read-depth=","normal-read-depth="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) < 1:
        usage()
        sys.exit(1)
    
    tumor_vcf_filenames  = []
    normal_vcf_filenames = []
    window_num   = None
    overlap_size = 2
    for opt, arg in opts:
        if opt == '-t':
            for filename in glob(f'{arg}/chr*.pileup.tsv'):
                tumor_vcf_filenames.append(filename)
        elif opt == '-n':
            for filename in glob(f'{arg}/chr*.pileup.tsv'):
                normal_vcf_filenames.append(filename)
        elif opt == '-w': window_num   = int(arg)
        elif opt == '-s': 
            print(arg)
            overlap_size = int(arg)
        elif opt in ['--tumor-read-depth']:  tumor_avg_depth = float(arg)
        elif opt in ['--normal-read-depth']: normal_avg_depth = float(arg)
    
    if len(tumor_vcf_filenames) == 0 and len(normal_vcf_filenames) == 0:
        print("VCF file error.")
        sys.exit(1)
    if len(normal_vcf_filenames) > 0 and len(tumor_vcf_filenames) > 0 and len(normal_vcf_filenames) != len(tumor_vcf_filenames):
        print("Matched normal VCF files have different number from Tumor VCF files.")
        sys.exit(1)

    if len(normal_vcf_filenames) > 0 and len(tumor_vcf_filenames) > 0:
        for tumor_vcf_filename in tumor_vcf_filenames:
            # get current chrName
            tmp = tumor_vcf_filename.split('/')[-1].split('.')
            chrom = ['chr' in i for i in tmp]
            chrName = [tmp[i] for i in range(len(chrom)) if chrom[i]][0]

            normal_vcf_filename = [s for s in normal_vcf_filenames if f'{chrName}.' in s][0]
            print(tumor_vcf_filename, normal_vcf_filename)
            # getCN&AF
            tumor_CN, tumor_AF, normal_CN, normal_AF = get_tumorCN_matchedNormal(tumor_vcf_filename, tumor_avg_depth, 
                normal_vcf_filename, normal_avg_depth, window_num=window_num, overlap_size=overlap_size)
            # write to individual files
            writeCN2file(tumor_CN, f'{chrName}.tumor.CN.tsv')
            writeCN2file(tumor_AF, f'{chrName}.tumor.AF.tsv')
            writeCN2file(normal_CN, f'{chrName}.normal.CN.tsv')
            writeCN2file(normal_AF, f'{chrName}.normal.AF.tsv')

    elif len(tumor_vcf_filenames) > 0:
        for filename in tumor_vcf_filenames:
            print(filename)
            # getCN&AF
            tumor_CN,tumor_AF = get_CN_data(filename, tumor_avg_depth, window_num=window_num, overlap_size=overlap_size)
            writeCN2file(tumor_CN, f'{chrName}.tumor.CN.tsv')
            writeCN2file(tumor_AF, f'{chrName}.tumor.AF.tsv')

    elif len(normal_vcf_filenames) > 0:
        for filename in normal_vcf_filenames:
            print(filename)
            # getCN&AF
            normal_CN,normal_AF = get_CN_data(filename, normal_avg_depth, window_num=window_num, overlap_size=overlap_size)
            writeCN2file(normal_CN, f'{chrName}.normal.CN.tsv')
            writeCN2file(normal_AF, f'{chrName}.normal.AF.tsv')
    else:
        print('Error with pileup VCF file directory input.')
        sys.exit(1)

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