def parse_ctgpaf(line, ctglen_thres, ctgspan_thres, mq_thres):
    '''
    Return a list of attributes from a line of paf file. Filter out if the total contig length,
    the aligned contig length or mapping quality is shorter than a certain threshold.
    ------
    Input: 
        line: one row of a given paf file
        ctglen_thres: the total contig length threshold, default is 50kb.
        ctgspan_thres: the aligned contig length threshold, default is 20kb.
        mq_thres: the mapping quality threshold, default is 1.
    Output:
        paf_line: a python list containing 10 selected attributes from line.
    '''
    tmp = line.rstrip().split('\t')
    indices = [*list(range(9)), 11]
    # 0ctgName, 1ctgLen, 2ctgStart, 3ctgEnd, 4ctgStrand, 5chrName, 6chrLen, 7chrStart, 8chrEnd, 9MQ
    paf_line = [int(tmp[_]) if tmp[_].isdigit() else tmp[_] for _ in indices]
    # if alignment is of low quality, filter out
    if paf_line[1] < ctglen_thres or (paf_line[3]-paf_line[2]+1)<ctgspan_thres or paf_line[9]<mq_thres or '_' in paf_line[5]:
        return None
    return paf_line

def get_hap(paf_filename):
    '''Retrieve hap1 or hap2 from paf file name.'''
    tmp = paf_filename.split('.')
    hap = ['hap' in i for i in paf_filename.split('.')]
    return [tmp[i] for i in range(len(hap)) if hap[i]][0]

def sort_ctg_in_chrom(ctg_pafs):
    '''sort a list of paf_line of contigs mapped to the same chromosome. 
    first sort by chrStart, then by ctgStart'''
    return sorted(ctg_pafs, key=lambda x: (x[7], x[2]))

def get_overlapRatio(A_start, A_end, B_start, B_end):
    '''calculate the overlap ratio between to sequence blocks A and B
    '''
    union   = max(A_end, B_end) - min(A_start, B_start) #useless variable...
    overlap = min(A_end, B_end) - max(A_start, B_start)
    if overlap <= 0:
        # There is NO overlap
        return 0
    return overlap / (A_end - A_start)

def ctg_not_in_censat(ctg_paf, chr_censat_dict, overlap_thres = 0.9):
    '''return True if contig alignment block is not in censat regions (default overlap>0.9), False otherwise.'''
    if len(chr_censat_dict)>1:
        chrName = ctg_paf[5]
        censat_regions = chr_censat_dict[chrName]
        A_start, A_end = ctg_paf[7], ctg_paf[8]

        for censat in censat_regions:
            B_start, B_end = censat
            B_start = min(B_start, abs(B_start-500000))
            B_end = B_end + 500000
            if get_overlapRatio(A_start, A_end, B_start, B_end) > overlap_thres:
                return False
    return True

def ctg_pass_multimap_filter(ctg_paf, ctg_mapcnt_dict, chr_censat_dict, len_thres = 25000):
    '''Return True if current contig alignment block is a valid candidate with translocation info,
    False otherwise. 
    Criteria: 
    1) spans across at least two chromosome; 2) not in centromeric regions; 
    3) the block length is at least 25000 bp by default, but ctglen_thres is always passed through
    4) if the alignment block is in telomere and subtelomeric regions
    5) if the alignment block is in acrocentric p-arm of chr13,14,15,21,22'''
    flag = False
    if ctg_paf[5] not in ctg_mapcnt_dict[ctg_paf[0]] and ctg_not_in_censat(ctg_paf, chr_censat_dict):
        if ctg_paf[3]-ctg_paf[2] > len_thres and ctg_paf[8] > 500000 and ctg_paf[8] < ctg_paf[6] - 500000:
            flag = True

    return flag

def split_by_chrom_helper(paf_filename, bed_filename, ctglen_thres, ctgspan_thres, mq_thres):
    '''helper to split paf file rows by chromosomes'''
    out_dict = defaultdict(list)
    chr_censat_dict = defaultdict(list)
    ctg_mapcnt_dict = defaultdict(list)
    ctg_pafs = []
    with open(paf_filename,'r') as paf_file:
        for line in paf_file.readlines():
            paf_line = parse_ctgpaf(line,ctglen_thres,ctgspan_thres,mq_thres)
            if paf_line:
                ctg_pafs.append(paf_line)
            else:
                continue

    if bed_filename:
        with open(bed_filename,'r') as bed_file:
            for line in bed_file.readlines():
                bed_line = line.rstrip().split('\t')[0:3]
                chr_censat_dict[bed_line[0]].append([int(bed_line[1]),int(bed_line[2])])

    for ctg_paf in ctg_pafs:
        out_dict[ctg_paf[5]].append(ctg_paf)
        # filter for ctg_mapcnt_dict
        if ctg_pass_multimap_filter(ctg_paf, ctg_mapcnt_dict, chr_censat_dict, len_thres = ctgspan_thres):
            ctg_mapcnt_dict[ctg_paf[0]].append(ctg_paf[5])
    
    rem_list = []
    for ctgName, chrList in ctg_mapcnt_dict.items():
        if len(chrList) < 2:
            rem_list.append(ctgName)
    [ctg_mapcnt_dict.pop(key) for key in rem_list]
    return out_dict, ctg_mapcnt_dict

def split_by_chrom(paf_filename, bed_filename=None, ctglen_thres = 50000, ctgspan_thres = 50000, mq_thres = 2):
    '''
    Input:
        paf_filename: paf file of contigs aligned to reference genome (e.g, CHM13 or hg38)
        bed_filename: bed file of difficult regions on the reference genome (e.g, cetromeric regions of CHM13)
        ctglen_thres: threshold that filters out contigs with short length
        ctgspan_thres: threshold that filters out contig alignment blocks with short length
        mq_thres: threshold that filters out contig alignment blocks with low mapping quality
    '''
    hap = get_hap(paf_filename)
    out_dict, ctg_mapcnt_dict = split_by_chrom_helper(paf_filename, bed_filename, ctglen_thres,ctgspan_thres,mq_thres)
    # write each one to files
    for chrom, ctg_pafs in out_dict.items():
        out_filename = chrom+'.'+hap+'.tsv'
        rows = sort_ctg_in_chrom(ctg_pafs)
        with open(out_filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            csvwriter.writerows(rows)
    
    # write contigs that might contain translocation info to one file
    with open(f'{hap}.ctg.multimap.tsv','w') as f:
        for ctgName, chrList in ctg_mapcnt_dict.items():
            for chrName in chrList:
                f.write(f'{ctgName}\t{chrName}\n')
    return None

def split_by_chrom_utg(paf_filename, bed_filename, utglen_thres = 20000, utgspan_thres = 10000, mq_thres = 2):
    out_dict, ctg_mapcnt_dict = split_by_chrom_helper(paf_filename, bed_filename, utglen_thres,utgspan_thres,mq_thres)
    # write each one to files
    for chrom, ctg_pafs in out_dict.items():
        out_filename = f'{chrom}.utg.tsv'
        rows = sort_ctg_in_chrom(ctg_pafs)
        with open(out_filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            csvwriter.writerows(rows)
   
    with open(f'utg.multimap.tsv','w') as f:
        for ctgName, chrList in ctg_mapcnt_dict.items():
            for chrName in chrList:
                f.write(f'{ctgName}\t{chrName}\n')
    return None

def usage():
    print("Usage: pafSplit.py -f paf_filename -b bed_filename")
    print("  -l INT      contig/unitig length threshold")
    print("  -s INT      contig/unitig span threshold")
    print("  -m INT      mapping quality threshold")
    print("  -f STR      paf file to be split")
    print("  -b STR      bed file containing difficult regions")
    print("Note: .")

def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:],"l:s:m:f:b:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) < 1:
        usage()
        sys.exit(1)

    # get paf_filenames, bed_filename, and parameters if specified
    ctgpaf_filename = utgpaf_filename = None
    bed_filename    = None
    len_thres = span_thres = mq_thres = 0
    l_flag = s_flag = m_flag = False
    for opt, arg in opts:
        if opt == '-f':
            if 'ctg' in arg:
                ctgpaf_filename = arg
            if 'utg' in arg:
                utgpaf_filename = arg
        elif opt == '-b':
            bed_filename = arg
        elif opt == '-l': 
            print(f'Length threshold is set to: {arg}')
            len_thres = int(arg)
            l_flag = True
        elif opt == '-s': 
            print(f'Span threshold is set to: {arg}')
            span_thres = int(arg)
            s_flag = True
        elif opt == '-m': 
            print(f'Mapping quality threshold is set to: {arg}')
            mq_thres = int(arg)
            m_flag = True
        
    print(f'Running {argv[0]} on contig file(s) {ctgpaf_filename}, unitig file {utgpaf_filename}')
    if ctgpaf_filename is not None:
        if not l_flag: len_thres = 50000
        if not s_flag: span_thres = 50000
        if not m_flag: mq_thres = 2
        split_by_chrom(ctgpaf_filename, bed_filename, ctglen_thres = len_thres, ctgspan_thres = span_thres, mq_thres = mq_thres)
    if utgpaf_filename is not None:
        if not l_flag: len_thres = 20000
        if not s_flag: span_thres = 10000
        if not m_flag: mq_thres = 2
        split_by_chrom_utg(utgpaf_filename, bed_filename, utglen_thres = len_thres, utgspan_thres = span_thres, mq_thres = mq_thres)

if __name__ == "__main__":
    import os
    import csv
    import sys
    import getopt
    from collections import defaultdict
    from collections import Counter
    main(sys.argv)
