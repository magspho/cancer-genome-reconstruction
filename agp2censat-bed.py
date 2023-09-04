# not useful
# providing a list of hg38 agp file from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP/
# extract the chromosome name from the file name and get centromeric positions:
# from the end of prev GL to the start of next GL, KI is centromeric regions

def extract_cenpos(agp_filename, idx_key = '# Format: '):
    idx = dict()
    with open(agp_filename,'r') as file:
        for line in file.readlines():
            if line.startswith('#'):
                if line.startswith(idx_key):
                    parsed_line = line.strip('#').strip().split()[1:]
                    if tsv_line is None:
                        continue
                    elif tsv_line[0] == 'CC' and tsv_line[1] == 'PR':
                        idx = dict()
                        for i in range(2,len(tsv_line)):
                            idx[tsv_line[i]] = i-2
                continue
            

def main(argv):
    #opts, args = getopt.getopt(argv[1:],"f:l:g:", 
    #                           ["paf=","len-thres=","gap-thres="])
    if len(argv[1:]) < 1:
        print("Usage: agp2censat-bed.py [filenames]")
        print("This simple script takes an agp file, extracts centromeric positions for each chromosome, and converts it to a simple three column bed file.")
        sys.exit(1)
    
    print(argv[1])
    for filename in glob(argv[1]):
        print(filename)
        chrName = filename.split('/')[-1].split('.')[0]
        
        out_filename = '.'.join(in_filename.split('.')[:-1])+'.fasta'
        with open(in_filename,'r') as file:
            with open(out_filename,'w') as fasta_file:
                for line in file.readlines():
                    tmp = line.rstrip().split('\t')
                    if len(tmp) != 2:
                        print(line)
                        print("The above line have more than 2 columns.")
                        sys.exit(1)

                    fasta_file.write(f">{tmp[0]}\n")
                    seq = tmp[1]
                    for i in range(0, len(seq), 60):
                        fasta_file.write(f"{seq[i:i+60]}\n")

if __name__ == "__main__":
    import sys
    
    main(sys.argv)