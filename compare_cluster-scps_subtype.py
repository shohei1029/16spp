#!/usr/bin/env python3

import sys

# Created on 2016.10.31 by Shohei Nagata.
# SCPSの結果ファイルを読み込ませ，クラスタリング結果とHIV-1もともとのサブタイプとの一致率を出す.
#  一致率は累計一致率的ななにかが最も高いやつを出す。クラスターとサブタイプとの対応が.
#  サブタイプ情報は，ID名に入ってるという前提の下でパースしてくる.

#若干hard coging部分ありにつき(サブタイプ名の数値化のところ)，subtype A,B,Cのみ対応

def get_cluster_info_and_subtype_from_scps_file(in_text):
    # {ID: [cluster number, subtype] }
    id_cl_st_d = {}
    for line in in_text:
        if '%' in line:
            continue
        if len(line) < 3: #空白行っぽいの飛ばす
            continue
        line = line.rstrip()

        cluster_id, cluster_num = line.split()
        st = cluster_id.split(':')[-1].split('-')[0]
        #integrade sub-subtypes
        if len(st) > 1:
            st = st[0]

        id_cl_st_d.update({cluster_id: [cluster_num, st]})

    return id_cl_st_d

def replace_by_table(subject:str, replace_table:dict):
    #とりま入力されてくる対象が一文字だと仮定する
    for k, v in replace_table.items():
        if subject == k:
            return str(v)

def calc_concordance_rate(id_cl_st_d):
    total_num = len(id_cl_st_d)
    concordance_rate_l = []

    #temp, hardcoding
    replace_tables = [ {'A':1, 'B':2, 'C':3}, {'A':2, 'B':3, 'C':1}, {'A':3, 'B':1, 'C':2} ] 

    for rep_table in replace_tables:
        counts = 0
        for v in id_cl_st_d.values():
            cl, st = v
            st = replace_by_table(st, rep_table)
    
            if cl == st:
                counts += 1
    
        rate = counts / total_num * 100
        concordance_rate_l.append(rate)

    return max(concordance_rate_l)


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as in_fh:
        id_cl_st_d = get_cluster_info_and_subtype_from_scps_file(in_fh)
    conc_rate = calc_concordance_rate(id_cl_st_d)
    print(conc_rate)
        
        
        






