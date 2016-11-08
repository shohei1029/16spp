#!/usr/bin/env python3

import os
import subprocess
import mimetypes #ファイルタイプ識別, テキスト，写真，etc.
import argparse
import time
import multiprocessing as mp


# 2016.5.10 by Shohei Nagata
# ディレクトリを指定し，その中のsimファイル全てに対してclusterxをかけてくれる．
# 分割数は2-10くらいで行う予定．
# ディレクトリ内一覧を取得した後，simかどうかの判別を行い(とりまsplitして3つだったら,って感じで)，simのみのリストを作る．
# 出力は，analysis内のclusterxフォルダに，inputと同じディレクトリ作成する．

# 2016.10.27
# 分割数の指定もオプションでできるように。"-"を使って，範囲指定もできる.

# 2016.11.08
# mpの使用。指定プロセス数のみのコマンドを投げる

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", type=str, help="Target Directory", required=True)
parser.add_argument("-c", "--num_clusters", type=str, default=3, help="sets the desired number of clusters to. you can use '-' to designate range of nums.")
args = parser.parse_args()

def find_all_files(directory):
    for root, dirs, files in os.walk(directory):
        #yield root
        for file in files:
            yield os.path.join(root, file)

def is_sim(file_path):
    if mimetypes.guess_type(file_path)[0] == 'text/plain':
        with open(file_path, 'r') as in_fh:
            for line in in_fh:
                l = line.split(' ')
                if len(l) == 3:
                    return True
                else:
                    return False
    else:
        return False

#def is_target(file_path, target_elems: '[[ドメイン達], [piリスト]]' ):
#hardcoding
def is_target(file_path):
    """
    末尾に'_'を付けてから，file_pathに含まれているかみる仕様に。
    """
#関数外に切り離したい。が，filter()使うときに関数に対象リスト以外の引数を与える方法がわからない。
    t_doms = "RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C".split(' ')
    #t_pis_num = list(range(65,85)) #2016.11.07
    t_pis_num = list(range(1,65)) #2016.11.08
    t_pis = ['pi'+str(x) for x in t_pis_num]
    target_elems = [t_doms, t_pis]
#

    file_name = file_path.split('/')[-1]
    flags = 0
    for elemlists in target_elems:
        for elem in elemlists:
            elem = str(elem) + '_'
            if elem in file_name:
                flags += 1
    if flags > len(target_elems):
        print("something wrong. check the target terms.")
    elif flags == len(target_elems):
        return True
    else:
        return False

def run_clusterx(input_file, output_dir, num_clusters):
    input_file_only = input_file.split('/')[-1]
    if '-' in str(num_clusters):
        cl_range_start = int(num_clusters.split('-')[0])
        cl_range_end = int(num_clusters.split('-')[1])
        for i in range(cl_range_start, cl_range_end):
            cmd = "nohup clusterx -c {c} -o {dir}/scps-c{c}_{fn} {inf} &".format(c=i,dir=output_dir,fn=input_file_only,inf=input_file)
            subprocess.call(cmd, shell=True)
            print(cmd)
            time.sleep(2)
    else:
        i = int(num_clusters)
        cmd = "nohup clusterx -c {c} -o {dir}/scps-c{c}_{fn} {inf} &".format(c=i,dir=output_dir,fn=input_file_only,inf=input_file)
        subprocess.call(cmd, shell=True)
        print(cmd)
        time.sleep(20)
#        print(cmd)

def run_cmd(cmd):
    subprocess.call(cmd, shell=True)

def run_clusterx_mp(input_file, output_dir, num_clusters, num_proc):
    input_file_only = input_file.split('/')[-1]
    cmd_l = []
    if '-' in str(num_clusters):
        cl_range_start = int(num_clusters.split('-')[0])
        cl_range_end = int(num_clusters.split('-')[1])
        for i in range(cl_range_start, cl_range_end):
            cmd = "clusterx -c {c} -o {dir}/scps-c{c}_{fn} {inf}".format(c=i,dir=output_dir,fn=input_file_only,inf=input_file)
            cmd_l.append(cmd)
    else:
        i = int(num_clusters)
        cmd = "clusterx -c {c} -o {dir}/scps-c{c}_{fn} {inf}".format(c=i,dir=output_dir,fn=input_file_only,inf=input_file)
        cmd_l.append(cmd)

    pool = mp.Pool(num_proc)
    outs = pool.map(run_cmd, cmd_l)
    return outs

    #funcをつかう

if __name__ == '__main__':
    files = list(find_all_files(args.input_dir))
    files_sim = list(filter(is_sim, files))

#暫定,hardcoding
    files_sim = list(filter(is_target, files_sim))

    output_dir = "../analysis/clusterx/{}".format(args.input_dir.rstrip('/').split('/')[-1])
    os.makedirs(output_dir, exist_ok=True)

    for simfile in files_sim:
#        run_clusterx(simfile, output_dir, args.num_clusters)

        #multiprocessing
        num_proc = mp.cpu_count()
        run_clusterx_mp(simfile, output_dir, args.num_clusters, num_proc)

