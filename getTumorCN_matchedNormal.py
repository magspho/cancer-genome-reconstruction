# get copy number of tumor from matched normal

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
    rd_gt =  sum([int(_) for _ in read_gt])
    h1_gt =  sum([int(_) for _ in hap1_gt])
    h2_gt =  sum([int(_) for _ in hap2_gt])
    read_depth = [int(_) for _ in parsed_line[col2idx['read']].split(':')[1].split(',')]
    
    if rd_gt == 2 or read_depth[0] <= depth_thres or read_depth[1] <= depth_thres:
    	return None
    # **Note**: if one cancer chromosome has odd number of copy number, 
    # h1_gt or h2_gt might be 1 since one haplotype from h2 may go to h1...?
    elif h1_gt == 1 or h2_gt == 1:
        return None
    elif h1_gt == h2_gt:
        print(f"Missing heterozygosity at {parsed_line[col2idx['#CHROM']]}: {parsed_line[col2idx['POS']]}")
        return None
    else:
        h1_depth = read_depth[int(h1_gt/2)]
        h2_depth = read_depth[int(h2_gt/2)]
        return [h1_depth, h2_depth]

def matchedNormalGenerator(normal_vcf_filename):
	'''Since normal_read_depth is ordered by variant positons, we call curr_nrd using this generator'''
	for curr_nrd in normal_read_depth:
		yield curr_nrd

def get_hap_depth_matchedNormal(parsed_line, curr_nrd, depth_thres=3):
	'''If we have a matched normal (mode_tumor=True), then we find the correct 
	hifiasm-assembled tumor haplotype with the help of corresponding heterozygous variant 
    from the matched normal even in case of hemizygosity in tumor sample.'''
	read_gt = parsed_line[col2idx['read']].split(':')[0].split('/')
    hap1_gt = parsed_line[col2idx['hap1']].split(':')[0].split('/')
    hap2_gt = parsed_line[col2idx['hap2']].split(':')[0].split('/')
    # gap in read
    if '.' in read_gt or ('.' in hap1_gt and '.' in hap2_gt):
        return None

    def process_gt(gt):
    	'''internal function to deal with gaps in genotypes'''
    	if '.' in gt: return 0
    	return sum([int(_) for _ in gt])
    rd_gt =  process_gt(read_gt)
    h1_gt =  process_gt(hap1_gt)
    h2_gt =  process_gt(hap2_gt)
    read_depth = [int(_) for _ in parsed_line[col2idx['read']].split(':')[1].split(',')]
    
    if rd_gt == 2 or read_depth[0] <= depth_thres or read_depth[1] <= depth_thres:
    	curr_pos = int(parsed_line[col2idx['POS']])

    # **Note**: if one cancer chromosome has odd number of copy number, 
    # h1_gt or h2_gt might be 1 since one haplotype from h2 may go to h1...?
    elif h1_gt == 1 or h2_gt == 1:
        return None
    elif h1_gt == h2_gt:
        print(f"Missing heterozygosity at {parsed_line[col2idx['#CHROM']]}: {parsed_line[col2idx['POS']]}")
        return None
    else:
        h1_depth = read_depth[int(h1_gt/2)]
        h2_depth = read_depth[int(h2_gt/2)]
        return [h1_depth, h2_depth]

def read_depth_from_vcf(vcf_filename, normal_vcf_filename = None, mode_tumor = False):
    '''
    Shared function between get_CN_data and get_tumorCN_matchedNormal with a mode_tumor switch.
	Note the variant could be heterozygous or hemizygous in tumor sample with a matched normal.
	Input:
		vcf_filename:		STR 	a pileup VCF file of one chromosome from the tumor sample
		normal_read_depth:	LIST	the position of a variant, hap1 read depth, and hap2 read depth
									Output of read_depth_from_vcf from the matching normal VCF file
		mode_tumor:			Boolean	Default False. If True, the argument must include normal_read_depth
	Output:
		a list of read depth (position, hap1_depth, hap2_depth)
    '''
    read_depth = []
    if normal_vcf_filename and mode_tumor:
    	curr_nrd = matchedNormalGenerator(normal_vcf_filename)
    else:
    	print("Error: mode_tumor with normal_read_depth.")
    	sys.exit(1)
    for line in open(vcf_filename, "r"):
        parsed_line = line.strip("\n").split("\t")
        if len(parsed_line) < 2:
            continue
    	elif not mode_tumor:
        	hap_depth = get_hap_depth(parsed_line)
        # mode_tumor
        else:
        	hap_depth = get_hap_depth_matchedNormal(parsed_line, curr_nrd)
        if hap_depth:
            read_depth.append([int(parsed_line[col2idx['POS']]), *hap_depth])
    return read_depth

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

def get_CN_data(vcf_filename, avg_depth, window_num = 250, overlap_size=2):
	'''
	Without a matched normal, CN is calculated from the position of heterozygous SNP, 
	hap1 read depth, and hap2 read depth from a pileup VCF file.
	Input:
		vcf_filename:	STR		a pileup VCF file of one chromosome from a sample (tumor/normal).
		avg_depth: 		INT 	the HiFi sequencing read depth divided by ploidy level.
		window_num:		INT 	number of windows to bin the CN across one chromosome. Default is 250.
		overlap_size:	INT 	size overlap for window sliding when calculating CN within a bin.
								Default is 2, meaning the window slides 1/2 of bin size.
	Output:
		CN: a 3D list containing the index of each bin, the average CN of Hap1 and Hap2 in each bin, respectively.
	'''
    read_depth = read_depth_from_vcf(vcf_filename)
    CN = smooth_copy_number(read_depth, window_num, avg_depth, overlap_size)
    CN = np.array(CN).T
    return list(CN)

def get_tumorCN_matchedNormal(tumor_vcf_filename, avg_depth, normal_vcf_filename, window_num = 250, overlap_size=2):
	'''
	Tumor sample CN calculated with a matched normal. 
	New Input:
		normal_read_depth: the output of read_depth_from_vcf from the matching normal VCF file
	'''
    # get current chrName
    tmp = tumor_vcf_filename.split('/')[-1].split('.')
    chrom = ['chr' in i for i in tmp]
    chrName = [tmp[i] for i in range(len(chrom)) if chrom[i]][0]

    read_depth = read_depth_from_vcf(tumor_vcf_filename, normal_vcf_filename, mode_tumor=True)
    CN = smooth_copy_number(read_depth, window_num, avg_depth, overlap_size)
    CN = np.array(CN).T
    
    return list(CN)

def usage():
	print("Usage: getTumorCN.py --tumor_read_depth [HiFi read depth of tumor sample] -t [directory to tumor pileup_byChr files]")
    print("Options:")
    print("  -t STR      Directory containing pileup VCF files for one chromosome of tumor sample")
    print("  -n STR      Directory containing pileup VCF files for one chromosome of matched normal sample")
    print("  --tumor-read-depth		INT      HiFi read depth of tumor sample")
    print("  --normal-read-depth	INT      HiFi read depth of matched normal sample")
    print("Note: processed pileup VCF file split by chromosomes should have the format <chrName>.pileup.tsv")
    print("Provide matched normal sample to resolve CNV in hemizygous regions.")

def main(argv):
	try:
    	opts, args = getopt.getopt(argv[1:],"t:n:", 
                               ["tumor-read-depth=","normal-read-depth="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) < 1:
        usage()
        sys.exit(1)
    
    tumor_vcf_filenames = []
    normal_vcf_filenames = []
    for opt, arg in opts:
        if opt == '-t':
            for filename in glob(f'{arg}/chr*.pileup.tsv'):
                tumor_vcf_filenames.append(filename)
        elif opt == '-n':
            for filename in glob(f'{arg}/chr*.pileup.tsv'):
                normal_vcf_filenames.append(filename)
        elif opt in ['--tumor-read-depth']:  tumor_read_depth = arg
        elif opt in ['--normal-read-depth']: normal_read_depth = arg
    
    if len(normal_vcf_filenames) > 0:
        for filename in normal_vcf_filenames:
            print(filename)
            # getCN


    if len(tumor_vcf_filenames) > 0:
        for filename in tumor_vcf_filenames:
            print(filename)
            # getCN
    else:
        print('Error with tumor pileup VCF file directory input.')
        sys.exit(1)

if __name__ == "__main__":
	import os
	import sys
	import csv
    import getopt
    from glob import glob
	from statistics import mean
	import matplotlib.pyplot as plt
	from collections import defaultdict
    main(sys.argv)