#!/usr/bin/env python3

import sys
import statistics

# 2016.4.15
# Dmean (配列間の距離の平均) を計算するスクリプト
# BLASTの総当たりbitscoresをもとにしたsim(x,y)を1-sim(x,y)して非類似度として使う。
# simファイルを読みこませる。
# simファイルは，重複とか，双方向からのエッジが書いてない前提！
# なんでもいいから何か引数があれが，Dmeanのみを表示する。他スクリプトとの連携がしやすいように。

def cal_Dmean(dists, num_isolates):
    S = num_isolates
    return sum(dists) / (S*(S+1)/2)

if __name__ == '__main__':
    dists = []
    ids = set()
    for line in sys.stdin:
        line = line.rstrip()
        id1, id2, sim = line.split(' ')
        dissim = 1 - float(sim)
        dists.append(dissim)
        ids.add(id1)
    
    num_isolates = len(ids)
    variance = statistics.variance(dists)
    dmean = cal_Dmean(dists, num_isolates)

    if sys.argv:
        sys.stdout.write(str(dmean))
    else:
        print("number of ids: ", num_isolates, file=sys.stdout)
        print("Variance of distance: ", variance, file=sys.stdout)
        print("Dmean: ", dmean, file=sys.stdout)
