
import re
from collections import defaultdict
import time
import os

from Bio import SeqIO
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt

from progressbar import ProgressBar

# 2016.2.9 
# Shohei Nagata
#
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
        
                def alter_domain_name_for_multihit(dn): #暫定：マルチヒットが3回までだと仮定
                    if dn in self.domain_dict[query_acc].keys():
                        print("multihit!",dn)
                        dn += "_"
                        if dn in self.domain_dict[query_acc].keys():
                            dn += "_"
                            return dn
                        else:
                            return dn
#                        alter_domain_name_for_multihit(dn) #ほんとは再帰で何回マルチヒットでも大丈夫にしたい
                    else:
                        return dn 

#                print(self.domain_dict[query_acc][domain_name])
                tmp_dn = alter_domain_name_for_multihit(domain_name)
                domain_name = tmp_dn

                self.domain_name_dict[domain_name] = domain_des
                self.domain_dict[query_acc][domain_name]["start"] = domain_start
                self.domain_dict[query_acc][domain_name]["end"]   = domain_end

    def read_fasta(self, fasta_file):
        with open(fasta_file, 'r') as fasta_fh:
            self.fasta_records = list(SeqIO.parse(fasta_fh, "fasta"))

    def prep_data(self):
        print("loading files..")
        self.get_domain_region_from_hmmscan_outfile(self.hmmscan_outfile)
        self.read_fasta(self.aligned_fasta)

    def create_image_matrix(self):
        print("creating image matrix..")
        self.num_matrix = []
        pbar = ProgressBar(max_value=len(self.fasta_records)) #for just indicate progress bar
        for i,seq in enumerate(self.fasta_records): #計算量...orz
            if print_details: print("\n",seq.id)
            row = []
            domain_region_pos = set()

            for d_name,pos_dic in self.domain_dict[seq.id].items(): #計算量... 下のfor文の中で一緒に処理できるようにしたひ
                gapped_start = convert_seqpos_to_gapped_seqpos(pos_dic["start"],seq.seq)
                gapped_end = convert_seqpos_to_gapped_seqpos(pos_dic["end"],seq.seq)
                if print_details: print(d_name,gapped_start,gapped_end)

                tmp_domain_region_pos = set(range(gapped_start-1,gapped_end))
                dom_pos_overlap = domain_region_pos.intersection(tmp_domain_region_pos)
                if print_details and dom_pos_overlap:
                    print("overlapping..",dom_pos_overlap)
                domain_region_pos = domain_region_pos.union(tmp_domain_region_pos)

            for p,c in enumerate(seq.seq): 
                if c == '-':
                    row.append(0) #gap
                elif p in domain_region_pos:
                    row.append(2) #domain region
                else:
                    row.append(1) #normal aa
            self.num_matrix.append(row)
           
            pbar.update(i)
        pbar.finish()


#        print(self.num_matrix)

    def draw_image_from_matrix(self, outfile, aspect=0.4, xaxis=False):
        print("drawing image..")
        # define the colormap
        cmap = plt.cm.jet
        cmaplist = ['white', 'gray', 'blue'] #問題点：例えば3カラー要素あるのに0,1しかないと，1なのに3番目のカラーになる
        #cmaplist = ['white', 'gray', 'blue', 'red']
        # create the new map
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

        plt.imshow( self.num_matrix, interpolation='none', aspect=aspect, cmap=cmap) 
        if xaxis:
            plt.tick_params(labelleft='off')
        else:
            plt.axis('off')
        plt.savefig(outfile, dpi = 800)

        
        
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
    print_details = False
    os.makedirs("../results/domain_pos",exist_ok=True)
    stime = time.time()
    hiv_pol = VisualiseProteinDomainRegion("../data/A,B,C_mafft-linsi_HIV-1-gM-noRs_pol-aa_v3.fasta", "../analysis/Pfam-hmmscan/PfamA29-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt")
#    hiv_pol = VisualiseProteinDomainRegion("../data/mafft-linsi_test.fasta", "../analysis/Pfam-hmmscan/PfamA29-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt")
    hiv_pol.create_image_matrix()
#    hiv_pol.draw_image_from_matrix("./test_HIV-1-gM-noRs_pol-aa_v3_aspect5.png",5) #test
#    hiv_pol.draw_image_from_matrix("../results/HIV-1-gM-noRs_pol-aa_v3_aspect0.1.png",0.1)
#    hiv_pol.draw_image_from_matrix("../results/HIV-1-gM-noRs_pol-aa_v3_aspect0.2.png",0.2)
#    hiv_pol.draw_image_from_matrix("../results/domain_pos/A,B,C_HIV-1-gM-noRs_pol-aa_v3_aspect0.3.png",0.3)
#    hiv_pol.draw_image_from_matrix("../results/domain_pos/A,B,C_HIV-1-gM-noRs_pol-aa_v3_aspect0.4_axis.png",0.4,True)
#    hiv_pol.draw_image_from_matrix("../results/domain_pos/A,B,C_HIV-1-gM-noRs_pol-aa_v3_aspect0.4.png",0.4,False)
#    hiv_pol.draw_image_from_matrix("../results/domain_pos/A,B,C_HIV-1-gM-noRs_pol-aa_v3_aspect0.5.png",0.5)
#    hiv_pol.draw_image_from_matrix("../results/domain_pos/A,B,C_HIV-1-gM-noRs_pol-aa_v3_aspect0.5.png",0.6)
#    hiv_pol.draw_image_from_matrix("../results/domain_pos/A,B,C_HIV-1-gM-noRs_pol-aa_v3_aspect0.8.png",0.8)
#    hiv_pol.draw_image_from_matrix("../results/HIV-1-gM-noRs_pol-aa_v3_aspect1.png",1)
    print(time.time() - stime, "[s]")

