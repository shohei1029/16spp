#!/usr/bin/env python3

import sys

from Bio import SeqIO

#2016.3.8
#HXB2のgag-pol配列を付けてアライメントしたfastaの，gap付き配列を読み込ませ，swissprotからのタンパク質領域とかの情報を入れ,
#アライメント後のfastaでのpositionを出力する。それだけ
#その後，collectDomain_byPfam_getpos_inAlignedFasta_createSecFasta.py のsectionごとに配列切り出す機能を使い，polのタンパク質ごとの配列を取得する。

#2016.3.20
#swiss protよりHXB2 gag-pol のドメイン(Pfamに基づく)をコピペし，ドメインごとの配列を取得する
# major update: ドメインごとの配列にしたfastaの作成まで行う

def seqpos_to_gapped_seqpos(pos,seq):
#    '''
#    ギャップなしFASTAファイルのpositionをギャップ有りFASTAファイルのpositionへ変換
#    convert char position in ungapped FASTA to gapped FASTA. (I'm sorry I can't use English well.)
#    '''
    c_count = 0
    gapped_pos = 0
    for c in seq:
        if c != '-':
            c_count += 1
        gapped_pos += 1
        if c_count == pos:
            return gapped_pos

def read_fasta_todict(fasta_file):
    with open(fasta_file,'r') as fa_fh:
        fa_dict = SeqIO.to_dict(SeqIO.parse(fa_fh, "fasta"))
    return fa_dict

def make_seqsec_fasta(seq_secs_ld, alnfasta_dict):
    fa2_dict = alnfasta_dict

    for k,v in fa2_dict.items(): 
        sys.stderr.write("sequence length with gap: {}\n".format(len(v.seq)))
        break

#strするときにpos-1で指定することを忘れずに！
    for seqsec_name, spos_l in seq_secs_ld.items():
        for k_acc in fa2_dict.keys():
            fasta_header = "{}-{}".format(k_acc, seqsec_name)
            sec_seq   = str(fa2_dict[k_acc].seq[spos_l[0]-1:spos_l[1]-1])
            sec_seq_nogap = sec_seq.replace('-','')
            outfa_fh.write(">{}\n{}\n".format(fasta_header,sec_seq_nogap))

def create_sec_pos_gpd(sec_ld_nogap, gpdseq):
    sec_ld_gap = {}
    for p, pl in sec_ld_nogap.items():
        s = seqpos_to_gapped_seqpos(pl[0],gpdseq)
        e = seqpos_to_gapped_seqpos(pl[1],gpdseq)
        sec_ld_gap[p] = [s,e]
    return sec_ld_gap


gpdseq = "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLREDLAFL-Q-GKARE--------F----------------S----------S---------------E-----------Q------T-----------------R------A-------------------N------------------------------S-------P---------------T-------R----------R------ELQVWG---RD---------N----N------S--P-------S-------E-----A-------G----AD----R----Q-------G-----T-----V-S----FN----FPQVTLWQRPLVTIKIGGQLKEALLDTGADDTVL-EEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPISPIETVPVKLKPGMDGPKVKQWPL-TEEKIKALVEICTEMEKEGKISKI-GPENPYNTPVFAIKKKD--STKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKL------NWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSK-----DLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWETWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGIIQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVLFLDGIDKAQDEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTGATVRAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKG-GIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVVIQ--DNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED-----------------"

#end_pos = len(gpdseq) #1744 #chainのときはIntegraseの末端をアライメント配列の末端へと広げた.ドメインのときはやらない
#sec_ld_nogap = {'Protease':[489, 587], 'p51 RT':[588, 1027], 'p15':[1028, 1147], 'Integrase':[1148, 1435]}
sec_ld_nogap = {'RVP':[493, 586], 'RVT_1':[650, 821], 'RVT_thumb':[828, 891], 'RVT_connect':[905, 1006], 'RNase_H':[1023, 1144], 'Integrase_Zn':[1155, 1192], 'rve':[1202, 1308], 'IN_DBD_C':[1368, 1415]}

if __name__ == '__main__':
    outfa_fh = sys.stdout
    seq_sec = create_sec_pos_gpd(sec_ld_nogap, gpdseq)
    fa_d = read_fasta_todict("../data/removed-HXB2_mafft-linsi_with-HXB2_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta")
    make_seqsec_fasta(seq_sec, fa_d)


