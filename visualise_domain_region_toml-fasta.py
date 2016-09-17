#!/usr/bin/env python3

import re
from collections import defaultdict
import time
import os

import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt

from Bio import SeqIO
import toml
from progressbar import ProgressBar

# 2016.2.9 
# Shohei Nagata
#
#アライメントされたfastaファイルと，hmmscanのdombtblout結果ファイルを読み込まs，予測された領域をdrawする。
#前提：hmmscanの結果ファイルに入っているモチーフ・ドメインの領域がちゃんとfastaファイル内の配列に存在すること
#検討：ドメイン別に色を変えるか

# 2016.5.11
# アライメントされたfastaファイルと，domain,region,chains等の位置が書かれたtomlファイルと読み込ませ，drawする. ように。
# tomlは，gap付きアライメントされたfastaを基に位置を変換されたやつを読みこませる。
#  変換スクリプトはconvert_toml-seqpos_toPosWithGap_1stSeqInFastaAsStandard.py

#アスペクトの自動調節ほしいなぁ


class VisualiseProteinSequenceRegion(object):
    def __init__(self, fasta_in, toml_file):
        self.aligned_fasta = fasta_in
        self.toml_file = toml_file
        self.prep_data()

    def read_fasta(self, fasta_file):
        with open(fasta_file, 'r') as fasta_fh:
            self.fasta_records = list(SeqIO.parse(fasta_fh, "fasta"))

    def read_toml(self, toml_file):
        with open(toml_file) as tomlfh:
            self.seqregion_dict = toml.loads(tomlfh.read())

    def prep_data(self):
        print("loading files..")
        self.read_fasta(self.aligned_fasta)
        self.read_toml(self.toml_file)

    def create_image_matrix(self):
        print("creating image matrix..")
        self.num_matrix = []
        pbar = ProgressBar(max_value=len(self.fasta_records)) #for just indicate progress bar

        seqregions_poss = set()
        for name, l in self.seqregion_dict.items():
            seqregions_poss.add(range(l[0]-1, l[1]-1))

        for i,seq in enumerate(self.fasta_records): #計算量...orz
            if print_details: print("\n",seq.id)
            row = []

            for p,c in enumerate(seq.seq): 
                if c == '-':
                    row.append(0) #gap
                elif p in seqregions_poss:
                    row.append(2) #domain region
                else:
                    row.append(1) #normal aa
            self.num_matrix.append(row)
           
            pbar.update(i)
        pbar.finish()

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
    os.makedirs("../results/region_pos",exist_ok=True)
    stime = time.time()

    aspect = 0.1
    test = VisualiseProteinSequenceRegion("../data/removed-HXB2_mafft-linsi_with-HXB2_HIV-1-gM-A,B,C-noRs_env-aa.fasta", "./toml/alignd_HIV-1_env_regions_HXB2.toml")
    test.create_image_matrix()
    test.draw_image_from_matrix("../results/region_pos/HIV-1-gM-A,B,C-noRs_env-aa_aspect{asp}_axis.png".format(asp=aspect),aspect,True)
    print(time.time() - stime, "[s]")

