# gfa2seq: 
# less -S K562-r2.hic.r_utg.gfa | grep S > /hlilab/yujieguo/K562_analysis/K562-r2.hic.r_utg.seq
# awk '/^S/{print ">"$2"\n"$3}' in.gfa | fold > out.fa

def main(argv):
    #opts, args = getopt.getopt(argv[1:],"f:l:g:", 
    #                           ["paf=","len-thres=","gap-thres="])
    if len(argv[1:]) < 1:
        print("Usage: seqidFile2fasta.py [seqid.tsv] [contigs.fa]")
        print("This simple script takes a seqid tsv file containing scaffold information and")
        print("contigs fasta file output from hifiasm and convert it to scaffold fasta file.")
        sys.exit(1)
    
    in_filename = argv[1]
    in_filename = in_filename.split('/')[-1]
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