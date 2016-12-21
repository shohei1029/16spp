
import os.path

import numpy as np
import pandas as pd
import networkx as nx

from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation

np.set_printoptions(threshold=np.inf)

# Created on 2016.12.21 by Shohei Nagata.
# simファイルとSCPS結果ファイルを読み込ませる
# 行列の形で持つようにし，scikit-learnによる解析(他クラスタリング手法，統計量算出等)を行えるようにする


def read_from_sim(in_sim_file):
    G = nx.Graph()
    with open(in_sim_file) as in_fh:
        for line in in_fh:
            line = line.rstrip()
            #source, edge_attr, target = line.split(" ")
            source, target, edge_attr = line.split(" ")
            G.add_edge(source, target, weight=float(edge_attr) )
    return G

def draw_network(G, node_color="r"):
    import matplotlib #to set use('Agg') 
    matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
    import matplotlib.pyplot as plt
#    import seaborn as sns
#    sns.set_style("whitegrid", {'grid.linestyle': '--'})
#    sns.set_palette("YlGnBu")

    nx.draw_networkx(G, with_labels=False, node_color=node_color, cmap=plt.cm.YlGnBu)
    plt.savefig('out_networkx.png')


if __name__ == "__main__":
    in_sim_file = "/Users/NagataShohei/Documents/bio-study/16spp/analysis/sim/A,B,C_pi92_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt"
    in_pkl = in_sim_file + ".pkl"
    if os.path.isfile(in_pkl):
        print("reading pickled sim data..")
        G = nx.read_gpickle(in_sim_file + ".pkl")
    else:
        read_from_sim(in_sim_file)
        nx.write_gpickle(G, in_sim_file + ".pkl")


    # G_adjmat = nx.to_numpy_matrix(G, nonedge=np.nan)
    #adjmat_df = nx.to_pandas_dataframe(G, nonedge=np.nan)
    adjmat_df = nx.to_pandas_dataframe(G)

#    clf = SpectralClustering(n_clusters=3, affinity="precomputed")
#    clf = KMeans(n_clusters=3)
    print("clustering..")
    clf = AffinityPropagation(affinity="precomputed")

    cluster_label = clf.fit_predict(adjmat_df)

    print("drawing..")
    draw_network(G, node_color=cluster_label)






#print(G_adjmat_df)

