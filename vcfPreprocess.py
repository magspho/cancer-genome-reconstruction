from io import StringIO
import sys
col2idx = None

def list_remove(inputList, indicesList):
    indicesList = sorted(indicesList, reverse=True)
    for indx in indicesList:
        # checking whether the corresponding iterator index is less than the list length
        if indx < len(inputList):
            # removing element by index using pop() function
            inputList.pop(indx)
        
def difficult_region(bed_file):
    '''
    A generator to iterate through a bed file containing GRCh38 difficult regions 
    which will be filtered out in the vcf file during preprocessing step.
    '''
    for line in open(bed_file,'r'):
        if line.startswith("#"):
            continue
        else:
            parsed_line = line.strip("\n").split("\t")
            # yield chrName and start and end of position range
            yield parsed_line[0], int(parsed_line[1]), int(parsed_line[2])

def check_dr(chrom, pos, prev_flt, dr_filter):
    '''Check if currect variant position falls in difficult regions'''
    if dr_filter is None:
        # Not checking (e.g. in the case of CHM13)
        return False, None
    flt = prev_flt
    is_in_dr = False
    while flt:
        dr_chrom, dr_start, dr_end = flt
        # determine whether dr_chrom is in front of or behind chrom
        if chrom != dr_chrom:
            dr_chrNum = dr_chrom.strip('chr')
            chrNum = chrom.strip('chr')
            if dr_chrNum.isnumeric() and chrNum.isnumeric():
                if int(dr_chrNum) < int(chrNum):
                    flt = next(dr_filter,None)
                    continue
                else: break
            elif dr_chrNum.isnumeric():
                flt = next(dr_filter,None)
                continue
            elif chrNum.isnumeric(): break
            else:
                if dr_chrNum < chrNum or (chrNum != 'X' and chrNum != 'Y'):
                    flt = next(dr_filter,None)
                    continue
                else: break
        # chrom == dr_chrom
        else:
            if dr_end <= pos:
                flt = next(dr_filter,None)
                continue
            elif dr_start > pos:
                break
            # dr_start <= pos and dr_end > pos
            else:
                is_in_dr = True
                break
    prev_flt = flt
    return is_in_dr, prev_flt

def get_readkey(col2idx, words):
    readkey = None
    for key, value in col2idx.items():
        flag = True
        for i in range(len(words)):
            if words[i] not in key:
                flag = False
        if flag:
            print(key,value)
            readkey = key 
    return readkey
            

def read_vcf(file_name, ref_bp = 1, alt_num = 1, qual_thres = 5, depth_thres = 10, bed_file=None, out_file = None):
    '''
    A function to read in vcf files from htsbox-pileup-r345 and filter out 
    genome position with low quality with specified threshold.
    This is for later to assess the phasing qualities of K562 assembled by hifiasm
    
    Input:
        - file_name: vcf file input name
        - ref_bp: The maximum length of bases of the reference allele.
                  Default is 1, indicating SNPs. If >1, then microsatellite is also considered.
        - alt_num: 
                  The maximum number of alternative alleles called. Default is 1.
        - qual_thres: 
                  Filtering threshold for the mapping quality score of ALT (cnt of read depth). Default is 5.
        - depth_thres:
                  The minimum reads depth sum at the reference position. Default is 30.
    '''
    global col2idx
    print(f"Running read_vcf({file_name})...")
    skipped_idx = []
    if bed_file:
        dr_filter = difficult_region(bed_file)
        prev_flt  = ['chr1',0,0]
    else:
        dr_filter = prev_flt = None
    line_count = pass1_count = pass2_count = pass3_count = pass4_count = 0
    readkey = None
    for line in open(file_name, "r"):
        if line.startswith("##"):
            continue
        parsed_line = line.strip("\n").split("\t")

        # Save col2idx -- won't be using ID, FILTER, INFO, FORMAT
        if line_count == 0 and line.startswith("#CHROM"):
            col2idx = {}
            lst = ['ID', 'FILTER', 'INFO', 'FORMAT']
            for i in range(len(parsed_line)):
                if parsed_line[i] in lst:
                    skipped_idx.append(i)
                    continue
            list_remove(parsed_line, skipped_idx)
            for i in range(len(parsed_line)):
                colName = parsed_line[i]
                if '/' in colName:
                    colName = colName.split('/')[-1].strip('.bam')
                if 'normal' in colName and 'read' in colName: # we have matched normal reads in the vcf in addition to tumor reads
                    readkey = colName
                col2idx[colName] = i
            # if no matched normal, col2idx has 8 entries; else 11 entries with the normal reads and haps.
            # col2idx = {'#CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3, 'QUAL': 4, 'read': 5, 'hap1': 6, 'hap2': 7}
            if readkey is None:
                readkey = get_readkey(col2idx,['read'])
            keys = col2idx.keys()
            print(keys)
            out_file.write('\t'.join(keys) + '\n')

        # Skip if encounter empty lines
        elif len(parsed_line) < 2:
            continue
            
        # Saving actual data
        else:
            line_count  += 1
            list_remove(parsed_line, skipped_idx)
            
            chrom = parsed_line[col2idx['#CHROM']]
            pos   = int(parsed_line[col2idx['POS']])
            ref   = parsed_line[col2idx['REF']]
            alt   = parsed_line[col2idx['ALT']]
            qual  = int(parsed_line[col2idx['QUAL']])
            
            rd_info  = parsed_line[col2idx[readkey]]
            depths   = [int(i) for i in rd_info.split(':')[1].split(',')]
            rd_depth = sum(depths)
            
            is_in_dr, prev_flt = check_dr(chrom, pos, prev_flt, dr_filter)
            can_write = True
            # pass if 1) #bp of REF > ref_bp or #num of ALT > alt_num or #bp of ALT > ref_bp (***scrutiny needed)
            if len(ref) > ref_bp or any(len(alt_snp) > ref_bp for alt_snp in alt.split(',')):
                # exclude len(alt.split(',')) > alt_num
                pass1_count += 1
                can_write = False
            # pass if 2) QUAL < qual_thres
            elif qual < qual_thres:
                pass2_count += 1
                can_write = False
            # pass if 3) Allelic depth for both the ref and alt alleles < 30
            elif rd_depth < depth_thres:
                pass3_count += 1
                can_write = False
            # pass if 4) the position belongs to difficult regions
            elif is_in_dr:
                pass4_count += 1
                can_write = False
            else:
                out_file.write('\t'.join(parsed_line) + '\n')
    print(f"Total number of variants: {line_count}.")
    print(f"Number of variants with more than {ref_bp} ref bases or {alt_num} alt alleles: {pass1_count}; {pass1_count/line_count*100:.3f}%.")
    print(f"Number of remaining variants with quality < {qual_thres}: {pass2_count}; {pass2_count/line_count*100:.3f}%.")
    print(f"Number of remaining variants with read depth < {depth_thres}: {pass3_count}; {pass3_count/line_count*100:.3f}%.")
    if bed_file:
        print(f"Number of remaining variants that are in difficult regions: {pass4_count}; {pass4_count/line_count*100:.3f}%.\n")

argv = sys.argv
if len(argv) < 2:
    print('Usage: vcfPreprocess.py [pileup_vcf] [difficult_region_bed]')
    sys.exit(1)
if len(argv) == 3:
    bed_file = argv[2]
vcf_file = argv[1]

out_file_name = f"processed.nodifficultregion.{vcf_file}"
with open(out_file_name, 'w') as f:
    read_vcf(vcf_file, bed_file=bed_file, out_file=f)