
import re
from collections import defaultdict
import time

from Bio import SeqIO
import matplotlib #to set use('Agg') 
import matplotlib.pyplot as plt
matplotlib.use('Agg')

class VisualiseProteinDomainRegion(object):
    def __init__(self, fasta_in, hmmscan_in):
        self.domain_dict = defaultdict(lambda:defaultdict(lambda:defaultdict()))
        self.domain_name_dict = defaultdict(str) #set的な用途で。
        self.aligned_fasta = fasta_in
        self.hmmscan_outfile= hmmscan_in

        self.prep_data()

    def get_domain_region_from_hmmscan_outfile(self, hmmscan_file):
        with open(hmmscan_file, 'r') as hmmscan_fh:
            for line in hmmscan_fh:
                if line[0] == '#':
                    continue
                line = line.rstrip()
        
                p = re.compile('\s+')
                hmmscandata = p.split(line)
        
                domain_name  = hmmscandata[0]
                #query_acc    = hmmscandata[3].split('|')[0]
                query_acc    = hmmscandata[3]
                domain_start = int(hmmscandata[17])
                domain_end   = int(hmmscandata[18])
                domain_des   = (' ').join(hmmscandata[22:])
        
                self.domain_name_dict[domain_name] = domain_des
                self.domain_dict[query_acc][domain_name]["start"] = domain_start
                self.domain_dict[query_acc][domain_name]["end"]   = domain_end

    def read_fasta(self, fasta_file):
        with open(fasta_file, 'r') as fasta_fh:
            self.fasta_records = list(SeqIO.parse(fasta_fh, "fasta"))

    def prep_data(self):
        self.get_domain_region_from_hmmscan_outfile(self.hmmscan_outfile)
        self.read_fasta(self.aligned_fasta)

    def create_image_matrix(self):
        self.num_matrix = []
        row = []
        for seq in self.fasta_records: #計算量...orz
            for c in seq: 
                if c == '-':
                    row.append(0)
                else:
                    row.append(1)
            self.num_matrix.append(row)

        plt.imshow( self.num_matrix, interpolation='none' ) )
        plt.savefig("./test.png")
        
        

        



if __name__ == '__main__':
    stime = time.time()
    test = VisualiseProteinDomainRegion("../data/mafft-linsi_HIV-1-gM-noRs_pol-aa_v3.fasta", "../data/Pfam-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt")
    test.create_image_matrix()
    print(time.time() - stime, "[s]")

