#!/usr/bin/env python3

import sys
import os
from collections import defaultdict
import itertools
import subprocess
import re

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO

import sn_utils

np.set_printoptions(threshold=np.nan)

# 2016.8.16 by Shohei Nagata.
# CRE (np.array)の要素位置は，1からスタートにしてある.そのままアライメントの位置と対応している.
# ->若干平均値とかに影響しちゃいそうな気がする.みてみたけどぱっとみ変わらなかった。

# 2016.9.6
# class化


# 作業メモ
# ../analysis/CRE/内にディレクトリを作成,fastaファイル名をもとに？
# fastaとgroup sizeが記載されたファイルを読みこませる
# groupの組み合わせを計算し，組み合わせだけのfastaファイルを生成する。
# それらfastaファイルに対してhmmbuildを実行
# .hmmファイルをパース，位置となんかをとってきてnp.array形式かなんかに格納
# plotする。arrayをplotするそんな関数作ってた気がする。
#
# ちゃんとクラス作ってやらないとなぁ...(遠い目) -> 対応.

def calc_kl_divergence(p, q):
    return stats.entropy(p, q) #default base: e

def cal_kld_2darray(p_nda,q_nda):
    kld_ll = []
    if len(p_nda) != len(q_nda):
        print("error: different sequence length..",)
        quit()

    for i in range(0, len(p_nda)):
        kld_ll.append( calc_kl_divergence(p_nda[i], q_nda[i]) )

    return np.array(kld_ll)

def normalize_array(array):
    mu = np.mean(array)
    std = np.std(array)
    z = (array - mu) / std
    return z

def smooth_array(array, window_size, min_valid_fraction=0.30):
    ma_array = np.ma.array(array)
    smoothed_divergence = sn_utils.running_average_masked(ma_array, window_size, min_valid_fraction=min_valid_fraction)
    return smoothed_divergence


class CalcPlotCRE(object):
    def __init__(self, fasta_path, gs_file_path, toml_path=None):
        self.in_fasta  = fasta_path
        self.gs_file   = gs_file_path
        self.toml_file = toml_path

        name_core = self.in_fasta.split('/')[-1].split('.')[0]
        self.workdir = "../analysis/CRE/{}".format(name_core)

        if self.toml_file:
            self.annopos_d = sn_utils.read_tomlfile(self.toml_file)

    def parse_gs_file_gen_gs(self, gsf):
        with open(gsf, 'r') as fh:
            for line in fh:
                for gs in line.rstrip().split(' '):
                    yield int(gs)
    
    def group_fasta(self, fa, gs_g):
        gfasta_d = defaultdict(dict)
        with open(fa, 'r') as fah:
            gs = gs_g.__next__()
            cgs = gs
            gfasta_d[gs] = []
            for i, record in enumerate(SeqIO.parse(fah, 'fasta')):
                gfasta_d[gs].append(record)
                if i == cgs-1:
                    try:
                        gs = gs_g.__next__()
                        cgs += gs
                        gfasta_d[gs] = []
                    except StopIteration:
                        break
        return gfasta_d
#        self.gfasta_d = gfasta_d
    
    def create_allcombinations_fasta_files(self, gfasta_d, workdir): #returns output fasta file path list
        outfiles_path_l = []
        combis = itertools.combinations(gfasta_d.keys(), len(gfasta_d.keys())-1) #-> [(a,b),(a,c),(b,c)] ,3つの要素の場合
    #    print(list(combis));quit()
    
        #グループ単体配列のfasta
        self.combis_singles = list(gfasta_d.keys())
        for k_gs, v_rcd_l in gfasta_d.items():
            outpath_single = "{}/{}.fasta".format(workdir, k_gs)
            outfiles_path_l.append(outpath_single)
            with open(outpath_single, 'w') as out_fah:
                for rcd in v_rcd_l:
                    SeqIO.write(rcd, out_fah, 'fasta')
    
        #グループ組み合わせのfasta
        self.combis_l = []
        for combi in combis:
            merged_rcd_l = []
            outpath_combi = "{}/{}.fasta".format(workdir, str(combi).replace("(","").replace(")","").replace(" ",""))
            outfiles_path_l.append(outpath_combi)
            self.combis_l.append(str(combi).replace("(","").replace(")","").replace(" ",""))
            for gs in combi:
                merged_rcd_l.extend(gfasta_d[gs])
            with open(outpath_combi, 'w') as out_fah:
                SeqIO.write(merged_rcd_l, out_fah, 'fasta')

        return outfiles_path_l
#        self.outfiles_path_l = outfiles_path_l

    def prep_data(self):
        os.makedirs(self.workdir, exist_ok=True)
        gs_g = self.parse_gs_file_gen_gs(self.gs_file)
        gfasta_d = self.group_fasta(self.in_fasta, gs_g)
        self.all_fasta_files_l = self.create_allcombinations_fasta_files(gfasta_d, self.workdir)
#        print(fafiles_l)
#        print(self.combis_l)
#        print(self.combis_singles)
    
    def run_hmmbuild_fasta_in_paths(self, fapath_l):
        outfiles_path = []
        for fapath in fapath_l:
    #        cmd = "echo {outf} {inf}".format(outf=fapath.replace(".fasta", ".hmm"), inf=fapath) #for test
            cmd = "hmmbuild {outf} {inf}".format(outf=fapath.replace(".fasta", ".hmm"), inf=fapath)
            subprocess.call(cmd, shell=True)
            outfiles_path.append(fapath.replace(".fasta", ".hmm"))
        return outfiles_path

    def run_hmmbuild(self):
        self.hmm_outfiles = self.run_hmmbuild_fasta_in_paths(self.all_fasta_files_l)
        
    def parse_hmm(self, hmm_file):
        with open(hmm_file, 'r') as fh:
            flag = False
            hmmll_l = []
            for line in fh:
                line = line.rstrip()
                if flag:
                    ll = re.split("\s+", line)
                    if len(ll) == 27:
                        hmmll_l.append(ll)
                elif flag == False and "HMM " in line:
                    flag = True
                if "//" in line:
                    break
        return hmmll_l
    
    def extract_aa_p(self, hmmll_l): #hmll_l -> hmm形式のデータのとこで各サイトの初行が1行づつ入っている
        seqlength = int(hmmll_l[-1][-5])
        ln_p_na = np.zeros((seqlength+1,20))  #ドメイン領域のポジションとか1からスタートだから，それに合わせる. 0番目を開ける->正規化時の平均値等に影響が.
    
        #2から22番目までがアミノ酸出現頻度(ln)，-5番目がポジション
        for hmml in hmmll_l:
    #        print(hmml)
    #        print(hmml[2:22])
            ln_p_na[int(hmml[-5]) ][0:20] = hmml[2:22] #
    
        p_na = np.exp(ln_p_na) #log(ln)を外す
        return p_na
    
    def calc_CRE_from_hmm(self):
        kld_l = []
        for s in self.combis_singles:
            hmm_m = self.parse_hmm("{}/{}.hmm".format(self.workdir, s))
            p_m = self.extract_aa_p(hmm_m)
            for cb in self.combis_l:
                if str(s) not in cb.split(","):
                    hmm_o = self.parse_hmm("{}/{}.hmm".format(self.workdir, cb))
                    p_o = self.extract_aa_p(hmm_o)
            kld_l.append(cal_kld_2darray(p_m, p_o))

        self.cre = sum(kld_l)
        self.cre_z = normalize_array(self.cre)

    def extract_high_pos(self, cre_z_thres=3.0): #書籍&元論文では3.0
        highposs = np.where(self.cre_z >= cre_z_thres)
        for name, pos_l in self.annopos_d.items():
            highposs_within = highposs[0][np.where( (pos_l[0] <= highposs[0]) & (highposs[0] <= pos_l[1]) )] 
            if len(highposs_within) != 0:
                print("CRE (Z-score) >= {} @ {}".format(cre_z_thres, name), highposs_within)

    def extract_high_pos_local(self, cre_z_thres=3.0): #書籍&元論文では3.0
        highposs = np.where(self.cre_z >= cre_z_thres)
        for name, pos_l in self.annopos_d.items():
            highposs_within = highposs[0][np.where( (pos_l[0] <= highposs[0]) & (highposs[0] <= pos_l[1]) )] -pos_l[0]+1
            if len(highposs_within) != 0:
                print("CRE (Z-score) >= {} @@ {}".format(cre_z_thres, name), highposs_within)
    
    def plot_nparray_with_annopos(self, array, annopos_d, outfile="./out_CRE_annopos.png", ylabel="Cumulative relative entropy", width=8):
        #plot setting
        sns.set('poster', 'whitegrid', 'dark', font_scale=1.2,
                rc={"lines.linewidth": 2, 'grid.linestyle': '--'})
        fig_width = width #standard:8, wide:50
        font_size = 20
        #
    
        fig, axs = plt.subplots(2,1, 
                                sharex=True, 
                                #figsize=(fig_width, fig_width*0.6),
                                figsize=(fig_width, 4.8), #height 固定
                                gridspec_kw={'height_ratios':[6, 1]})
    
        #plot entropy
        ax = axs[0]
        ax.plot(array, c='k')
        ax.set_ylabel(ylabel, fontsize=font_size)
    
        ax = axs[1]
        sns.set_style('white')
        sn_utils.draw_annopos(ax, annopos_d)
        ax.set_yticks([])
        ax.set_xlabel('Position', fontsize=font_size)
    
        # Final touches
        plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.1, h_pad=0.5, w_pad=0.4)
        fig.savefig(outfile, dpi=180)
        sns.plt.close()

    def plot(self):
        self.plot_nparray_with_annopos(self.cre_z, self.annopos_d, outfile="{wd}/out_CREz_annopos.png".format(wd=self.workdir), ylabel="CRE (Z-score)", width=8)

    def plot_wide(self):
        self.plot_nparray_with_annopos(self.cre_z, self.annopos_d, outfile="{wd}/out_CREz_annopos_wide.png".format(wd=self.workdir), ylabel="CRE (Z-score)", width=50)
        


if __name__ == "__main__":
    #plot settings

    argvs = sys.argv
    in_fasta = argvs[1]
    gs_file = argvs[2]
    toml_file = argvs[3]

    test = CalcPlotCRE(in_fasta, gs_file, toml_file)
    #test = CalcPlotCRE(in_fasta, gs_file, "/Users/NagataShohei/Documents/bio-study/16spp/scripts/toml/HIV-1_gag-pol_Pfam_HXB2_+p6-pol.toml")
    test.prep_data()
#    test.run_hmmbuild()
    test.calc_CRE_from_hmm()
    test.extract_high_pos(3.0)
    test.extract_high_pos(1.5)
    test.extract_high_pos(1.0)
    test.extract_high_pos_local(3.0)
    test.extract_high_pos_local(1.5)
    test.extract_high_pos_local(1.0)
#    test.plot()
#    test.plot_wide()

#    ws = 50
#    plot_nparray_with_annopos(cre, annopos_d, outfile="{wd}/out_CRE_annopos.png".format(wd=workdir))
#    plot_nparray_with_annopos(cre_z, annopos_d, outfile="{wd}/out_CREz_annopos2.png".format(wd=workdir), ylabel="CRE (Z-score)", width=8)
#    plot_nparray_with_annopos(cre_z, annopos_d, outfile="{wd}/out_CREz_annopos2_wide.png".format(wd=workdir), ylabel="CRE (Z-score)", width=50)

#    for ws in 10,20,50,100,120,150,200:
#        cre_ma = smooth_array(cre_z, ws)
#        plot_nparray_with_annopos(cre_ma, annopos_d, outfile="{wd}/out_CREz_annopos_mv{w}_wide2.png".format(wd=workdir, w=ws))
