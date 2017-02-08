
import os.path

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AffinityPropagation, SpectralClustering, AgglomerativeClustering, DBSCAN 

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

def get_cluster_label_from_scps_file(in_scps_file):
    # {ID: cluster label }
    id_cl_d = {}
    with open(in_scps_file) as in_fh:
        for line in in_fh:
            if '%' in line:
                continue
            if len(line) < 3: #空白行っぽいの飛ばす
                continue
            line = line.rstrip()
    
            cluster_id, cluster_lab = line.split()
    
            id_cl_d.update({cluster_id: int(cluster_lab)})

    return id_cl_d

def sort_accto_list(kv_d, key_order_l):
    out_val_l = []
    #keyが全部dictに存在している前提. ひとまず.
    for k in key_order_l:
        out_val_l.append(kv_d[k])

    return out_val_l

def draw_network(G, node_size=100, node_color="r", outfile="out_networkx.png"):
#    import seaborn as sns
#    sns.set_style("whitegrid", {'grid.linestyle': '--'})
#    sns.set_palette("YlGnBu")

    print("drawing..")
    nx.draw_networkx(G, with_labels=False, node_size=node_size, node_color=node_color, cmap=plt.cm.Set1)
    plt.savefig(outfile, dip=180)
    plt.close()
    print("save network plot as {} ..".format(outfile))

def perform_clusterings_and_draw(adj_mat, n_clusters=3, outdir="./"):
    print("clustering.. and silhouette analyis..")

    #K-means -> shape=(n_samples, n_features)なため不適?
#    clf = KMeans(n_clusters=n_clusters)
#    cluster_label = clf.fit_predict(adjmat_df)
#    draw_network(G, outfile= outdir + "out_networkx_KMeans.png", node_color=cluster_label)

#    clf = AffinityPropagation(affinity="precomputed", verbose=True)
#    cluster_label = clf.fit_predict(adjmat_df)
#    draw_network(G, outfile= outdir + "out_networkx_AffinityPropagation.png", node_color=cluster_label)
#    silhouette_analysis(adj_mat, cluster_label, outfile= outdir + "out_silhouette_analysis_AffinityPropagation.png")

    clf = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', assign_labels='discretize')
    cluster_label = clf.fit_predict(adjmat_df)
    draw_network(G, outfile= outdir + "out_networkx_SpectralClustering-discretize_c{}.png".format(n_clusters), node_color=cluster_label)
    silhouette_analysis(adj_mat, cluster_label, outfile= outdir + "out_silhouette_analysis_SpectralClustering-discretize_c{}.png".format(n_clusters))

#    clf = AgglomerativeClustering(n_clusters=n_clusters, affinity="precomputed", linkage="average")
#    cluster_label = clf.fit_predict(adjmat_df)
#    draw_network(G, outfile= outdir + "out_networkx_AgglomerativeClustering-average.png", node_color=cluster_label)
#    silhouette_analysis(adj_mat, cluster_label, outfile= outdir + "out_silhouette_analysis_AgglomerativeClustering-average.png")
#
#    clf = AgglomerativeClustering(n_clusters=n_clusters, affinity="precomputed", linkage="complete")
#    cluster_label = clf.fit_predict(adjmat_df)
#    draw_network(G, outfile= outdir + "out_networkx_AgglomerativeClustering-complete.png", node_color=cluster_label)
#    silhouette_analysis(adj_mat, cluster_label, outfile= outdir + "out_silhouette_analysis_AgglomerativeClustering-complete.png")
#
#    clf = DBSCAN(metric="precomputed", n_jobs=4)
#    cluster_label = clf.fit_predict(adjmat_df)
#    draw_network(G, outfile= outdir + "out_networkx_DBSCAN.png", node_color=cluster_label)
#    #silhouette_analysis(adj_mat, cluster_label, outfile= outdir + "out_silhouette_analysis_DBSCAN.png")
#    # -> 1 clusterしかなくてエラー

def silhouette_analysis(X_adj_mat, y_cluster_label, outfile="./out_silhouette_analysis.png"):
    from sklearn.metrics import silhouette_samples

    print("silhouette analysis..") 
    cluster_labels = np.unique(y_cluster_label)
    n_clusters = cluster_labels.shape[0]
    # シルエット係数を計算
    silhouette_vals = silhouette_samples(X_adj_mat, y_cluster_label, metric='precomputed')
    y_ax_lower, y_ax_upper = 0, 0
    yticks = []
    for i, c in enumerate(cluster_labels):
        c_silhouette_vals = silhouette_vals[y_cluster_label == c]
        c_silhouette_vals.sort()
        y_ax_upper += len(c_silhouette_vals)
        color = plt.cm.jet(i / n_clusters)
        plt.barh(range(y_ax_lower, y_ax_upper),
                 c_silhouette_vals,
                 height=1.0,
                 edgecolor='none',
                 color=color)
        yticks.append((y_ax_lower + y_ax_upper) / 2) # クラスラベルの表示位置を追加
        y_ax_lower += len(c_silhouette_vals) # 底辺の値に棒の幅を追加

    silhouette_avg = np.mean(silhouette_vals)
    plt.axvline(silhouette_avg, color="red", linestyle="--") # 係数の平均値に破線をひく
#    plt.yticks(yticks, cluster_labels)
    plt.yticks(yticks, cluster_labels + 1)
    plt.ylabel("Cluster")
    plt.xlabel("Silhouette coefficient")
#    plt.show()
    plt.savefig(outfile, dip=180)
    plt.close()


if __name__ == "__main__":
    outdir = "../results/sklearn_clustering/A,B,C_pi92_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3/"
    in_sim_file = "/Users/NagataShohei/Documents/bio-study/16spp/analysis/sim/A,B,C_pi92_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt"
    in_pkl = in_sim_file + ".pkl"
    if os.path.isfile(in_pkl):
        print("reading pickled sim data..")
        G = nx.read_gpickle(in_sim_file + ".pkl")
    else:
        G = read_from_sim(in_sim_file)
        nx.write_gpickle(G, in_sim_file + ".pkl")

    # SCPS
    in_scps_file = "/Users/NagataShohei/Documents/bio-study/16spp/analysis/clusterx/HIV-1-gM-A,B,C-noRs_pol-aa/scps-c3_A,B,C_pi92_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt"
    cl_id_lab_d = get_cluster_label_from_scps_file(in_scps_file)
    scps_label = np.array( sort_accto_list(cl_id_lab_d, G.nodes()) )

    # Creating affinity matrix
    # G_adjmat = nx.to_numpy_matrix(G, nonedge=np.nan)
    #adjmat_df = nx.to_pandas_dataframe(G, nonedge=np.nan) #本来はエッジなしもあるからこっち。ただnanだとクラスタリング時にエラーが出てしまう。。
    adjmat_df = nx.to_pandas_dataframe(G)

     # 対称行列か確認 -> True
#    arr = adjmat_df.as_matrix()
#    print((arr.transpose() == arr).all())

    #SCPS
#    draw_network(G, outfile="out_networkx_scps.png", node_color=scps_label)
#    silhouette_analysis(adjmat_df, scps_label -1 , outfile=outdir + "out_silhouette_analysis_SCPS.png") #scikit-learnでのクラスターラベルは0から始まるようになってて合わせるため-1をする

    for i in range(2,11):
        perform_clusterings_and_draw(adjmat_df, n_clusters=i, outdir=outdir)

    










