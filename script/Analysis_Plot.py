#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 23:01:41 2021

@author: mac
"""

# Analysis-Plot1
def K_Choice_KNNDensity(PCAData,k=20):
    N = PCAData.shape[0]
    starttime = datetime.datetime.now()
    Distance = np.triu(np.round(squareform(pdist(np.array(PCAData), metric='euclidean')),4),1)
    Distance = Distance.flatten()
    Distance = Distance[Distance>0]
    neighbor_dist = np.round(np.quantile(Distance,0.01),4)
    
    X = []
    Y = []
    for k in [10,20,30,40,50,60,70,80,90,100]:
        KNN_Dist = kneighbors_graph(PCAData, n_neighbors=k, mode='distance', metric='minkowski', p=2, include_self=False)
        KNN_Neighbor_Ind  = np.array([np.array(KNN_Dist[i].nonzero()[1]) for i in range(N)])
        KNN_Neighbor_Dist = np.array([KNN_Dist[i, KNN_Neighbor_Ind[i]].toarray()[0] for i in range(N)])
        reciprocal  = [list(1/i) for i in KNN_Neighbor_Dist]
        Den = np.array([sum(i) for i in reciprocal])
        Den = Den/Den.sum()*100000
        Freq = np.histogram(Den,bins=50)[0]
        Prob = np.round(Freq/Freq.sum(),4)
        SI = Shannon_Index(Prob)
        plt.figure(figsize=(8,6),dpi=1200)
        plt.hist(Den,bins=50)
        plt.ylabel("Frequency")
        plt.xlabel("Density")
        plt.title("Histogram of KNN Density")
        plt.savefig(name + "/Hist_KNN_Den_" + str(k)+ ".png")
        plt.close("all")  
        X.append(k)
        Y.append(SI)
        
    plt.figure(figsize=(8,6),dpi=1200)
    plt.plot(X, Y, marker='o', mec='r', mfc='w')
    plt.ylabel("Shannon Index")
    plt.xlabel("K")
    plt.ylim(0,6)
    plt.title("Choice of K for KNN Density")
    plt.savefig(name + "/K_Choice_KNNDensity.png")
    plt.close("all")  