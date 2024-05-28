# This version tries to assemble chromosomes with translocations.

def quick_display(ctg_segs):
    '''just get rid of contained ctgs and print out the primary contigs'''
    prev = None
    output = []
    for ctg in ctg_segs:
        if iscontained(prev, ctg, threshold=0.05):
            continue
        else:
            output.append(ctg)
            prev = ctg
    return output

def ctgline2ID(ctg):
    '''write attributes of a contig to SeqID'''
    ctg_len        = f'ctg_len={ctg[1]}'
    ctg_range      = f'range={ctg[2]}-{ctg[3]}'
    ctg_strand     = f'strand={ctg[4]}'
    chr_annotation = f'{ctg[5]}:{ctg[7]}-{ctg[8]}'
    return ' | '.join([ctg[0], chr_annotation, ctg_range, ctg_strand, ctg_len])

def iscontained(ctg_seg1, ctg_seg2, threshold=0.01):
    '''
    check if ctg_seg2 is contained in ctg_seg1, the default threshold is 0.01
    ***Note: threshold means if less that 1% of the ctg_seg2 base is out of the seg1 range at either end,
    it will be counted as contained.
    '''
    # first ctg_seg
    if ctg_seg1 is None:
        return False
    # ctg aligned to different chromosomes or prev chrStart is larger
    if ctg_seg1[5] != ctg_seg2[5] or ctg_seg1[7] > ctg_seg2[7]:
        return False

    seg2_len = ctg_seg2[1]
    if ctg_seg1[7] <= ctg_seg2[7] and ctg_seg1[8]-ctg_seg2[8] > -seg2_len*threshold:
        return True
    elif ctg_seg1[8] >= ctg_seg2[8] and ctg_seg2[7]-ctg_seg1[7] > -seg2_len*threshold:
        return True
    else:
        return False

def check_ctg_sorted(ctg_segs):
    '''return 1 if ctg_segs list is sorted in ascending order by chrStart, 0 otherwise'''
    flag = 0
    i = 1
    while i < len(ctg_segs):
        if(ctg_segs[i][7] > ctg_segs[i - 1][7]):
            flag = 1
        i += 1
    return flag

def check_ctg_len(ctg, prev, len_thres=50000):
    '''return False if current contig is shorter than len_thres or is contained in the previous one.
    Note: short/contained contig will be saved to HapIDs as it may have important genotype info.'''
    # first ctg
    if prev is None: return True
    # short contig
    elif ctg[1] < len_thres: return False
    # contained contig
#     elif iscontained(prev, ctg, threshold=0.01): 
#         print('contained contig:')
#         print('\t'.join([str(i) for i in prev]))
#         print('\t'.join([str(i) for i in ctg]))
#         return False
    else: return True

def get_ctg_range(ctg_segs):
    first = last = ctg_block = None
    if ctg_segs[0][4] != ctg_segs[-1][4]:
        print(f'{ctg_segs[0][0]}: Strand not consistent between the first and last segments.')
    # get first ctg_seg (smallest ctgstart)
    first = ctg_segs[0]
    # get last ctg_seg (largest ctgend)
    ctg_segs = sorted(ctg_segs, key=lambda l:l[3],reverse=True)
    last = ctg_segs[0]
    # get the most prominent ctg_seg strand
    ctg_segs = sorted(ctg_segs, key=lambda l:l[3]-l[2],reverse=True)
    strand = ctg_segs[0][4]
    # get smallest chrStart
    ctg_segs_chrS = sorted(ctg_segs, key=lambda l:l[7])
    chrS = ctg_segs_chrS[0][7]
    # get lasrgest chrEnd
    ctg_segs_chrE = sorted(ctg_segs, key=lambda l:l[8],reverse=True)
    chrE = ctg_segs_chrE[0][8]
    
    ctgb_len = last[3]-first[2]
    ctgb_aln = chrE-chrS
    ctgb_len_ratio = round((last[3]-first[2])/first[1],3)
    ctgb_aln_ratio = round((chrE-chrS)/first[6],3)
    ctg_block = [first[0],first[1],first[2],last[3],strand,first[5],first[6],chrS,chrE,
                 ctgb_len,ctgb_aln,ctgb_len_ratio,ctgb_aln_ratio]
    return ctg_block

def connect_ctg_segs(ctg_segs, gap_thres = 50000):
    '''Stitching the contigs with the same name together within the given list of contig segments. 
    Mark the ones with strand difference and large gap in alignment (which won't be connected).'''
    # sort by the starting position of contig coordinates
    ctg_segs = sorted(ctg_segs, key=lambda l:l[2])
    output = []
    first = last = None
    
    # single contig:
    if len(ctg_segs) == 1:
        output.append(get_ctg_range(ctg_segs))
    
    # check gaps between contig segments and chromosome coordinates
    # segments without large gaps form into blocks
    else:
        block_idx = 0
        for i in range(len(ctg_segs)-1):
            prev = ctg_segs[i]
            curr = ctg_segs[i+1]
            if curr[2]-prev[3] > gap_thres or abs(min((curr[7]-prev[8]),(curr[8]-prev[7]))) > gap_thres:
                print(f'{curr[0]}{curr[4]}: Gap from {prev[3]} to {curr[2]} on {curr[5]}:{prev[8]}-{curr[7]}')
                output.append(get_ctg_range(ctg_segs[block_idx:i+1]))
                block_idx = i+1
        output.append(get_ctg_range(ctg_segs[block_idx:]))
    
    # check unmapped regions of the contig
    last=ctg_segs[-1]
    if last[1]-last[3] > gap_thres:
        if last[4] == '+':
            print(f'{last[0]}{last[4]}: Gap from {last[3]} to end on {last[5]}:{last[8]}-unknown')
        else:
            print(f'{last[0]}{last[4]}: Gap from {last[3]} to end on {last[5]}:unknown-{last[7]}')
    return output

def get_scaffold(ctg_byChr_filename, gap_thres=50000, len_thres=50000):
    '''Get list of contigs (name+description) that form a scaffold. 
    Each list element will be the ID in FASTA file.'''
    prev = is_longctg = None
    # the following three variables are all dictionary of lists
    SeqIDs_tmp = defaultdict(list)
    SeqIDs = dict()
    HapIDs_tmp = defaultdict(list)
    
    with open(ctg_byChr_filename,'r') as ctg_file:
        for line in ctg_file.readlines():
            tmp = line.rstrip().split('\t')
            ctg = [int(_) if _.isdigit() else _ for _ in tmp]
            is_longctg = check_ctg_len(ctg, prev, len_thres) # aka pass len_thres and not contained
            prev = ctg
            if is_longctg:
                SeqIDs_tmp[ctg[0]].append(ctg)
            else:
                HapIDs_tmp[ctg[0]].append(ctg)
    for ctg_name, ctg_segs in SeqIDs_tmp.items():
        SeqIDs[ctg_name] = connect_ctg_segs(ctg_segs, gap_thres)

    return SeqIDs, HapIDs_tmp, SeqIDs_tmp

def lprint(lst):
    '''A nice way to print list of 1D/2D'''
    if isinstance(lst[0],str):
        print('\t'.join([str(val) for val in lst]))
    else:
        for val in lst:
            print('\t'.join([str(i) for i in val]))
    return None

def dprint(dicto):
    '''A nice way to print dictionary with values of list'''
    for keys,values in dicto.items():
        if isinstance(values[0],str):
            print('\t'.join([str(val) for val in values]))
        else:
            for val in values:
                print('\t'.join([str(i) for i in val]))
    return None

def write2tsv_helper(ctg_segs, reversed=False, threshold=0.01, len_thres = 50000):
    primary_tsv_lst = []
    nonpmry_tsv_lst = []
    ctg_segs = sorted(ctg_segs, key=lambda l:l[7])
    prev = None
    for ctg in ctg_segs:
        if iscontained(prev, ctg, threshold=threshold):
            nonpmry_tsv_lst.append(['CO',*ctg])
        elif ctg[3]-ctg[2] < len_thres:
            continue
        else:
            primary_tsv_lst.append(['PR',*ctg])
            nonpmry_tsv_lst.append(['//'])
            nonpmry_tsv_lst.append(['PR',*ctg])
            prev = ctg
    if reversed:
        new_nonpmry_tsv_lst = []
        prev_idx = len(nonpmry_tsv_lst)
        primary_tsv_lst.reverse()
        for line in primary_tsv_lst:
            idx = nonpmry_tsv_lst.index(line)
            new_nonpmry_tsv_lst.append(nonpmry_tsv_lst[idx])
            if prev_idx-1 != idx:
                new_nonpmry_tsv_lst.extend(nonpmry_tsv_lst[idx+1:prev_idx])
            prev_idx = idx
        nonpmry_tsv_lst = new_nonpmry_tsv_lst

    return primary_tsv_lst, nonpmry_tsv_lst

def write2tsv(SeqIDs, filedir, chrom, hap, no_trans=True, threshold=0.01, len_thres = 50000):
    '''
    write to two tsv files: 
    1. only primary contigs
    2. full file contained contigs as well in terms of reference genome coordinates
    '''
    file1_name = f'{filedir}/{chrom}.{hap}.primary.seqid.tsv'
    file2_name = f'{filedir}/{chrom}.{hap}.seqid.tsv'
    os.makedirs(filedir, exist_ok=True)
    file1 = open(file1_name, 'w', newline='')
    file2 = open(file2_name, 'w', newline='')
    # write headers of the full file
    ctg_header = '\t'.join(['ctgName','ctgLen','ctgStart','ctgEnd','ctgStrand',
                            'chrName','chrLen','chrStart','chrEnd',
                            'ctgbLen','chrbLen','ctgbRatio','chrbRatio'])
    header = '\n'.join(['CC\tDR\tchr\tstart\tend',
                       'CC\tCM\tp/q-arm',
                       f'CC\tPR\t{ctg_header}',
                       f'CC\tCO\t{ctg_header}\n'])
    file1.write(header)
    file2.write(header)
    
    if isinstance(SeqIDs, dict):
        ctg_segs = list(chain.from_iterable(SeqIDs.values()))
    else:
        ctg_segs = SeqIDs

    csvwriter1 = csv.writer(file1, delimiter='\t')
    csvwriter2 = csv.writer(file2, delimiter='\t')
    # only sort whole list if there's no translocation
    if no_trans:
        primary_tsv_lst, nonpmry_tsv_lst = write2tsv_helper(ctg_segs, threshold=threshold, len_thres = len_thres)
        csvwriter1.writerows(primary_tsv_lst)
        csvwriter2.writerows(nonpmry_tsv_lst)
    # if there's translocation, write it for each chromosome section
    else:
        sections = defaultdict(list)
        prev = ctg_segs[0]
        prev_idx = 0
        for idx in range(len(ctg_segs)):
            ctg = ctg_segs[idx]
            if prev[5] != ctg[5]:
                sections[prev[5]] = ctg_segs[prev_idx:idx]
                prev_idx = idx
            prev = ctg
        sections[prev[5]] = ctg_segs[prev_idx:]
        for chrName, ctg_segs in sections.items():
            if check_ctg_sorted(ctg_segs) == 1:
                primary_tsv_lst, nonpmry_tsv_lst = write2tsv_helper(ctg_segs)
            # if reverse sorted
            else:
                primary_tsv_lst, nonpmry_tsv_lst = write2tsv_helper(ctg_segs, reversed=True)
            csvwriter1.writerows(primary_tsv_lst)
            csvwriter2.writerows(nonpmry_tsv_lst)

    file1.close()
    file2.close()
    return None

# begins section of find translocations
def parse_tsvfile(filename):
    '''The tsvfile starts with 'CC', which indicates the legends for each row.
    A row that starts with 'PR' stands for primary contig,
    and a row that starts with 'CO' stands for contained contig.'''
    global idx
    ctg_segs = []
    with open(filename,'r') as tsv_file:
        for line in tsv_file.readlines():
            tsv_line = line.rstrip().split('\t')
            if tsv_line is None:
                continue
            elif tsv_line[0] == 'CC' and tsv_line[1] == 'PR':
                idx = dict()
                for i in range(2,len(tsv_line)):
                    idx[tsv_line[i]] = i-2
            elif tsv_line[0]=='PR' or tsv_line[0]=='CO':
                tsv_line = [int(_) if _.isdigit() else _ for _ in tsv_line]
                ctg_segs.append(tsv_line[1:])
            else:
                continue
    return ctg_segs

def stitch_translocation(filename_list, t_ctgName):
    print(filename_list)
    print(t_ctgName)
    ctg_segs_byChr = dict()
    t_idx_byChr = dict()
    out_ctg_segs = []
    for filename in filename_list:
        ctg_segs = parse_tsvfile(filename)
        curr_chrName = ctg_segs[0][idx['chrName']]
        ctg_segs_byChr[curr_chrName] = ctg_segs
    
    for chrom, ctg_segs in ctg_segs_byChr.items():
        for i in range(len(ctg_segs)):
            ctg = ctg_segs[i]
            if t_ctgName in ctg:
                t_idx_byChr[chrom]=i
    t_chroms = list(t_idx_byChr.keys())
    for j in range(len(t_chroms)-1):
        prev_ctg_segs = ctg_segs_byChr[t_chroms[j]]
        curr_ctg_segs = ctg_segs_byChr[t_chroms[j+1]]
        prev_t_idx = t_idx_byChr[t_chroms[j]]
        prev_t_ctg = prev_ctg_segs[prev_t_idx]
        curr_t_idx = t_idx_byChr[t_chroms[j+1]]
        curr_t_ctg = curr_ctg_segs[curr_t_idx]
        
        # this only works if one contig is only connecting one translocation
        if prev_t_ctg[idx['ctgStrand']] == curr_t_ctg[idx['ctgStrand']]:
            out_ctg_segs.extend(prev_ctg_segs[0:prev_t_idx+1])
            out_ctg_segs.extend(curr_ctg_segs[curr_t_idx:])
        else:
            out_ctg_segs.extend(prev_ctg_segs[0:prev_t_idx+1])
            out_ctg_segs.extend(reversed(curr_ctg_segs[0:curr_t_idx+1]))
    return out_ctg_segs

def stitch_translocation_byPos(filename_list, prev_pos, next_pos):
    '''
    pos: a tuple of chrName and chrPos
    '''
    return None

def get_filedir(filename):
    return '/'.join(filename.split('/')[0:-1])

def find_translocation(ctg_multimap_filename):
    '''
    After the output from get_scaffold is written to two tsv files for each normal chromosome, 
    we extract the contig names from the hap.ctg.multimap.tsv that indicates possible translocations.
    Then we find the corresponding filename_list that contain this ctgName. If not found in primary.seqid.tsv, 
    then seqid.tsv with contained contig blocks are considered.
    '''
    filedir = get_filedir(ctg_multimap_filename)
    ctg_mapcnt_dict = defaultdict(list)
    with open(ctg_multimap_filename,'r') as f:
        for line in f.readlines():
            tmp = line.rstrip().split('\t')
            ctg_mapcnt_dict[tmp[0]].append(tmp[1])
    for ctgName, chrList in ctg_mapcnt_dict.items():
        filename_list = []
        t_chr_arg = []
        if 'h1' in ctgName:
            hap = '1'
        elif 'h2' in ctgName:
            hap = '2'
        else:
            print(f"Error: check contig {ctgName}")
        for chrName in chrList:
            # stitch translocation from the seqid.tsv with contained contigs
            filename_list.append(f'{filedir}/assembly/{chrName}.hap{hap}.seqid.tsv')
            t_chr_arg.append(chrName.strip('chr'))
        ctg_segs = stitch_translocation(filename_list, ctgName)
        t_chr_arg = f'h{hap}_'.join(t_chr_arg)
        t_chr_arg = t_chr_arg+'h'+hap
        write2tsv(ctg_segs, filedir+'/assembly', 'chrt', t_chr_arg, no_trans = False)

def write2agp(SeqIDs, out_filename):
    '''subject to change'''
    with open(out_filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        # csvwriter.writerows(SeqIDs.values())
        for val in SeqIDs.values():
            csvwriter.writerows(val)
    return None

def write2fasta(fagz_file, SeqIDs, result_file):
    ctg_names = SeqIDs.keys()
    with gzip.open(fagz_file, "rt") as handle:
        with open(result_file, "w") as f:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in ctg_names:
                    tmp = SeqIDs[record.id]
                    # check if the contigs in seqids have different segments
                    if isinstance(tmp[0],str):
                        record.description = ctgline2ID(tmp)
                        SeqIO.write([record], f, "fasta")
                    else:
                        for seqid in tmp:
                            seqid_f = ctgline2ID(seqid)
                            record.description = seqid_f
                            SeqIO.write([record], f, "fasta")
            
def main(argv):
    opts, args = getopt.getopt(argv[1:],"d:l:g:", 
                               ["paf-dir=","len-thres=","gap-thres="])
    if len(opts) < 1:
        print("Usage: ctgSegment.py -d [directory to files]")
        print("Options:")
        print("  -d STR      Directory containing PAF files aligned between contigs and reference genome for one chromosome haplotype")
        print("  -l INT      len_threshold for including a contig alignment block in the scaffold. Default is 50k")
        print("  -g INT      gap_threshold for stitching two alignment block of the same contig together. Default is 50k")
        print("Note: processed ctg_aln2ref PAF file split by chromosome + haplotypes should have the format <chrName>.<haplotype>.tsv")
        sys.exit(1)
    
    paf_filenames = []
    len_thres = gap_thres = 50000
    for opt, arg in opts:
        if opt in ['-d','--paf-dir']:
            for filename in glob(f'{arg}/chr*.hap*.tsv'):
                paf_filenames.append(filename)
        elif opt in ['-l','--len-thres']: len_thres = int(arg)
        elif opt in ['-g','--gap-thres']: gap_thres = int(arg)
    
    # normal chromosomes
    filedir = get_filedir(paf_filenames[0])
    if len(paf_filenames) > 0:
        for filename in paf_filenames:
            print(filename)
            # extract chromsome and haplotype info
            tmp = filename.split('/')[-1].split('.')
            chrom = tmp[0]
            hap = tmp[1]
            seqid, hapid, tmp = get_scaffold(filename, len_thres=len_thres, gap_thres=gap_thres)
            write2tsv(seqid, filedir+'/assembly', chrom, hap)
    else:
        print('Error with PAF file input')
        sys.exit(1)
        
    # chromosomes with possible translocations
    ctg_multimap_filenames = [f'{filedir}/hap1.ctg.multimap.tsv',f'{filedir}/hap2.ctg.multimap.tsv']
    for filename in ctg_multimap_filenames:
        find_translocation(filename)

if __name__ == "__main__":
    import sys
    import os
    import getopt
    from glob import glob
    import csv
    import gzip
    import pandas as pd
    import numpy as np
    from Bio import SeqIO
    from itertools import chain
    from collections import defaultdict
    
    main(sys.argv)