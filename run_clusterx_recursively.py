#!/usr/bin/env python3

import os
import subprocess
import mimetypes #ファイルタイプ識別, テキスト，写真，etc.
import argparse
import time

# 2016.5.10 by Shohei Nagata
# ディレクトリを指定し，その中のsimファイル全てに対してclusterxをかけてくれる．
# 分割数は2-10くらいで行う予定．
# ディレクトリ内一覧を取得した後，simかどうかの判別を行い(とりまsplitして3つだったら,って感じで)，simのみのリストを作る．
# 出力は，analysis内のclusterxフォルダに，inputと同じディレクトリ作成する．

# 2016.10.27
# 分割数の指定もオプションでできるように。"-"を使って，範囲指定もできる.

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
        
def run_clusterx(input_file, output_dir, num_clusters):
    input_file_only = input_file.split('/')[-1]
    if '-' in num_clusters:
        cl_range_start = int(num_clusters.split('-')[0])
        cl_range_end = int(num_clusters.split('-')[1])
        for i in range(cl_range_start, cl_range_end):
            cmd = "nohup clusterx -c {c} -o {dir}/scps-c{c}_{fn} {inf} &".format(c=i,dir=output_dir,fn=input_file_only,inf=input_file)
            subprocess.call(cmd, shell=True)
            time.sleep(2)
    else:
        i = int(num_clusters)
        cmd = "nohup clusterx -c {c} -o {dir}/scps-c{c}_{fn} {inf} &".format(c=i,dir=output_dir,fn=input_file_only,inf=input_file)
        subprocess.call(cmd, shell=True)
#        print(cmd)


if __name__ == '__main__':
    files = list(find_all_files(args.input_dir))
    files_sim = list(filter(is_sim, files))

    output_dir = "../analysis/clusterx/{}".format(args.input_dir.rstrip('/').split('/')[-1])
    os.makedirs(output_dir, exist_ok=True)

    for simfile in files_sim:
        run_clusterx(simfile, output_dir, args.num_clusters)
