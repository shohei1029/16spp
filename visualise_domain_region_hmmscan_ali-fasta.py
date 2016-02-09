
import re
from collections import defaultdict
import time

from Bio import SeqIO
import matplotlib #to set use('Agg') 
import matplotlib.pyplot as plt
matplotlib.use('Agg')

from progressbar import ProgressBar

#アライメントされたfastaファイルと，hmmscanのdombtblout結果ファイルを読み込まs，予測された領域をdrawする。
#前提：hmmscanの結果ファイルに入っているモチーフ・ドメインの領域がちゃんとfastaファイル内の配列に存在すること

#検討：ドメイン別に色を変えるか

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
        self.pbar = ProgressBar(max_value=len(self.fasta_records)) #for just indicate progress bar
        for i,seq in enumerate(self.fasta_records): #計算量...orz
            row = []
            domain_region_pos = set()
            for d_name,pos_dic in self.domain_dict[seq.id].items(): #計算量... 下のfor文の中で一緒に処理できるようにしたひ
                gapped_start = convert_seqpos_to_gapped_seqpos(pos_dic["start"],seq.seq)
                gapped_end = convert_seqpos_to_gapped_seqpos(pos_dic["end"],seq.seq)
#                print(d_name,gapped_start,gapped_end)

                tmp_domain_region_pos = set(range(gapped_start-1,gapped_end))
#                print("overlapping..",domain_region_pos.intersection(tmp_domain_region_pos))
                domain_region_pos = domain_region_pos.union(tmp_domain_region_pos)

            for p,c in enumerate(seq.seq): 
                if c == '-':
                    row.append(0) #gap
                elif p in domain_region_pos:
                    row.append(2) #domain region
                else:
                    row.append(1) #normal aa
            self.num_matrix.append(row)
           
            self.pbar.update(i)


#        print(self.num_matrix)

    def draw_image_from_matrix(self, outfile):
        # define the colormap
        cmap = plt.cm.jet
        cmaplist = ['white', 'gray', 'blue'] #問題点：例えば3カラー要素あるのに0,1しかないと，1なのに3番目のカラーになる
        #cmaplist = ['white', 'gray', 'blue', 'red']
        # create the new map
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

        plt.imshow( self.num_matrix, interpolation='none', aspect=5, cmap=cmap) 
        plt.axis('off')
        plt.savefig(outfile, dpi = 500)

        self.pbar.finish()
        
        
def convert_seqpos_to_gapped_seqpos(pos,seq):
    '''
    ギャップなしFASTAファイルのpositionをギャップ有りFASTAファイルのpositionへ変換
    convert char position in ungapped FASTA to gapped FASTA. 
    '''
    c_count = 0
    gapped_pos = 0
    for c in seq:
        if c != '-':
            c_count += 1
        gapped_pos += 1
        if c_count == pos:
            return gapped_pos
        



if __name__ == '__main__':
    stime = time.time()
#    hiv_pol = VisualiseProteinDomainRegion("../data/mafft-linsi_HIV-1-gM-noRs_pol-aa_v3.fasta", "../data/Pfam-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt")
    hiv_pol = VisualiseProteinDomainRegion("../data/mafft-linsi_test.fasta", "../data/Pfam-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt")
    hiv_pol.create_image_matrix()
    hiv_pol.draw_image_from_matrix("./HIV-1-gM-noRs_pol-aa_v3.png")
    print(time.time() - stime, "[s]")

