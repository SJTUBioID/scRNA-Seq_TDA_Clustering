# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 21:29:16 2019

@author: shaolab
"""

## Import library

import os
import gc
import pandas as pd
import numpy as np
import matplotlib
import random
import seaborn as sns
import networkx as nx
import datetime
import itertools
#import dcor
import warnings


from math import pi
from copy import deepcopy
from time import time
from collections import Counter
from itertools import groupby
from scipy import stats
from scipy import sparse
from scipy import spatial as sp
from scipy.spatial import distance
from scipy.cluster import hierarchy as hc
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests

from sklearn.neighbors import NearestNeighbors
from sklearn import manifold
from sklearn.manifold import TSNE
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.manifold import Isomap
from sklearn.manifold import MDS
from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift
from sklearn.neighbors import KernelDensity
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import kneighbors_graph
from sklearn.model_selection import GridSearchCV
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.metrics import euclidean_distances
from sklearn.metrics import adjusted_mutual_info_score 
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import completeness_score
from sklearn.metrics import fowlkes_mallows_score
from sklearn.metrics import homogeneity_completeness_v_measure
from sklearn.metrics import homogeneity_score
from sklearn.metrics import mutual_info_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples

#from sklearn.metrics import jaccard_score
from sklearn.preprocessing import scale
from sklearn.preprocessing import MinMaxScaler
from umap import UMAP
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import colors as color
from matplotlib import gridspec
from matplotlib import cm
from matplotlib import ticker
from matplotlib.pylab import contour
from numpy.random import multivariate_normal
from matplotlib import pyplot as plt
from mpl_toolkits import axisartist
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors import kneighbors_graph

from scipy import sparse 
from scipy.special import gamma
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist



## Options
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 200)
warnings.filterwarnings("ignore")
matplotlib.use("Agg")
plt.style.use(['science','nature','no-latex'])
plt.rcParams['font.sans-serif'] = "Helvetica"



## color map
startcolor = '#FFFFFF' 
midcolor = 	'#008B8B'         
endcolor = '#8B0000'  
mycmap = color.LinearSegmentedColormap.from_list('MyCMAP',[startcolor,midcolor,endcolor])
cm.register_cmap(cmap=mycmap)


def Read10XGenomics(Path):
    if not Path.endswith("/"):
        Path = Path + "/"
    barcode = pd.read_csv(Path+"barcodes.tsv", header=None, index_col=None, sep="\t")
    
    if "genes.tsv" in os.listdir(Path):
        gene = pd.read_csv(Path+"genes.tsv", header=None, index_col=None, sep="\t")
    if "features.tsv" in os.listdir(Path):
        gene = pd.read_csv(Path+"features.tsv", header=None, index_col=None, sep="\t")
        
    mat = pd.read_csv(Path+"matrix.mtx", header=2, index_col=None, sep="\s+")
    barcode["Cell"] = barcode.index + 1
    barcode.columns = ["Barcode","Cell"]
    gene["Gene"] = gene.index + 1
    columns = list(gene.columns)
    columns[1] = "Symbol"
    gene.columns = columns
    mat.columns = ["Gene","Cell","UMI"]
    
    mat = pd.merge(gene,mat,on="Gene")
    mat = pd.merge(barcode,mat,on="Cell")
    mat = mat.loc[:,["Symbol","Barcode","UMI"]]
    
    mat.loc[:,["Symbol","Barcode"]] =  mat.loc[:,["Symbol","Barcode"]].astype('category')
    mat.loc[:,["UMI"]] = mat.loc[:,["UMI"]].astype("uint16")
    
    result = pd.pivot_table(mat,index="Symbol",columns="Barcode",values="UMI",aggfunc = np.sum, fill_value=0)
    result = result.fillna(0)
    return(result)

def UMI_Distribution_Plot(data,name):
    TotalUMI = data.sum(axis=0)
    plt.hist(TotalUMI,bins=50,color="coral",edgecolor="k")
    plt.xlabel("UMI")
    plt.title("UMI Distribution")
    plt.savefig(name + "/" + name + "_UMI-Distribution.png",dpi=300)
    plt.close("all")



##  Quality Control
def Quality_Control(data,cellNum=3,UMINumLowLimit=1000,UMINumUpLimit=10000,mtPct=0.1):
    # Select Genes
    data = data.loc[(data>0).sum(axis=1) >= cellNum,:]
    # Select Cells
    #data = data.loc[:,(data>0).sum(axis=0) >= geneNum]
    data = data.loc[:,data.sum(axis=0) >= UMINumLowLimit]
    #data = data.loc[:,data.sum(axis=0) <= UMINumUpLimit]
    # Select mitogenes
    mt_row = [r for r in data.index.tolist() if r.startswith("mt")]
    mt_percent=data.loc[mt_row].sum()/data.sum()
    data = data.loc[:, mt_percent <= mtPct]
    return(data)
    
def HighQuality_Data(data, UMINumLowLimit = 5000, cellNum=3, mtPct=0.2):
    # Select Cells
    data = data.loc[:,data.sum(axis=0) >= UMINumLowLimit]
    
    # DownSampling
    UMI = data.sum(axis = 0)
    for i in range(data.shape[1]):
        Prob = data.iloc[:,i]/UMI[i] 
        data.iloc[:,i] = np.random.multinomial(UMINumLowLimit, Prob)
        
    # Select Genes
    data = data.loc[(data>0).sum(axis=1) >= cellNum,:]    
    
    # Select mitogenes
    mt_row = [r for r in data.index.tolist() if r.startswith("mt")]
    mt_percent=data.loc[mt_row].sum()/data.sum()
    data = data.loc[:, mt_percent <= mtPct]        
    return(data)
#    
def Normalize_Rawdata(data):
    data = np.log(10000*(data/data.sum(axis=0))+1)
    #data = np.log2(100000*(data/data.sum(axis=0))+1)  #(SC3 + scTDA)
    return(data)


def Select_HEGenes(data,ngenes=500):
    data['std'] = data.std(axis=1)
    data = data.sort_values(by="std" , ascending=False)
    data = data.drop('std',axis=1,inplace=False)
    data = data.iloc[0:ngenes,]
    return(data)

def Select_HVGenes(data,ngenes=500):
    data['std'] = data.std(axis=1)
    data['mean'] = data.mean(axis=1)
    data['cv'] = data['std'] /data['mean']
    data = data.sort_values(by="cv" , ascending=False)
    data = data.drop('std',axis=1,inplace=False)
    data = data.drop('mean',axis=1,inplace=False)
    data = data.drop('cv',axis=1,inplace=False)
    data = data.iloc[0:ngenes,]
    return(data)

def Select_HEGL(data,ngenes=500):
    data['std'] = data.std(axis=1)
    data = data.sort_values(by="std" , ascending=False)
    data = data.drop('std',axis=1,inplace=False)
    data = data.iloc[0:500,]   
    
    Cormat = data.T.corr().abs()
    genelist = [Cormat.index.to_list()[0]]
    for i in range(1, Cormat.shape[0]):
        print(i)
        temp = Cormat.loc[genelist,Cormat.columns[i]]

        if temp.max() <= 0.1:
            genelist.append(Cormat.index[i])
        
    Result = data.loc[genelist,]
    return(Result)

def Scaling_Data(data):
    scaled = pd.DataFrame(scale(data,axis=1,with_mean=True,with_std=True),columns=data.columns,index=data.index)
    return(scaled)


def PCA_Permutation(data, name, nperm = 100, nPC = 50):
    Dimension = 100    
    Dimension = min(Dimension,data.shape[0])
    pca = PCA(n_components=Dimension, random_state = 1)
    PCA_s = pd.DataFrame(pca.fit_transform(data), index=data.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
    RealEigen = np.round(pca.explained_variance_,2)
    EigenPerm =  np.zeros((nperm,Dimension))
    
    for i in range(nperm):
        DataPerm = np.array(data.iloc[np.random.randint(low=0,high=data.shape[0],size=data.shape[0]),:])
        for j in range(DataPerm.shape[1]):
            np.random.shuffle(DataPerm[:,j])
        PCA_perm = pca.fit_transform(DataPerm)
        EigenPerm[i,] = pca.explained_variance_

    Eigen = pd.DataFrame(EigenPerm, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]        
    plt.figure(figsize=(8,3),dpi=600)
    x = np.linspace(1,Dimension,num=Dimension)
    y1 = RealEigen
    bar1 = plt.bar(x, y1, facecolor='coral',edgecolor='black', label = "EigenValue")
    for i in range(Eigen.shape[1]):
        y = Eigen.values[:,i]
        x = np.random.normal(i+1, 0.1, size=len(y))
        plt.plot(x, y, 'o', markerfacecolor='crimson', markeredgecolor='black',markeredgewidth=0.1, markersize=1.2)
    plt.xlabel("Principal Components")
    plt.xticks(np.linspace(0,nPC,num=nPC+1))
    plt.xlim((0,nPC))
    plt.ylabel("Eigen Values")
    plt.title(name)
    plt.savefig(name + "/PC_SelectionBarplot.pdf",dpi=600)
    plt.savefig(name + "/PC_SelectionBarplot.tiff",dpi=600)
    plt.close("all")
    Diff = list((RealEigen - Eigen.values.max(axis=0)) > 0)
    if False in Diff:
        N = Diff.index(False)
    else:
        N = Dimension
    return N


def tSNEPlot(tSNE,name,Label,legend={"loc":'center right', "bbox_to_anchor":(1.3,0.5),"fontsize":5, "ncol":2}):
    unique_labels = list(set(Label))
    unique_labels.sort()
    all_labels = deepcopy(unique_labels)
    
    if -1 in unique_labels:
        unique_labels.remove(-1)
    
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    plt.figure(figsize=(4,3))
    # Plot unassigned cells
    n = -1
    col = [0, 0, 0, 1]
    class_mask = [i == n for i in Label]
    tSNE_mask = tSNE.loc[class_mask,:]
    if tSNE_mask.shape[0] > 0:
        plt.plot(tSNE_mask.iloc[:,0],tSNE_mask.iloc[:,1],'o', markerfacecolor=tuple(col), markeredgecolor='k',markeredgewidth=0.2,markersize=2,label = str(n))

    # Plot assigned cells
    for n, col in zip(unique_labels, colors):
        #print(n,col)
        class_mask = [i == n for i in Label]
        tSNE_mask = tSNE.loc[class_mask,:]
        plt.plot(tSNE_mask.iloc[:,0],tSNE_mask.iloc[:,1],'o', markerfacecolor=tuple(col), markeredgecolor='k',markeredgewidth=0.2,markersize=2, label = str(n))
    

    plt.legend(labels = all_labels,loc=legend["loc"], bbox_to_anchor=legend["bbox_to_anchor"],fontsize=legend["fontsize"], ncol=legend["ncol"])
    N = len(unique_labels)
    #plt.title("tSNE, nCluster = " + str(N))
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    plt.subplots_adjust(right = 0.8) 
    #plt.savefig(name + "_tSNE.tiff",dpi=600)
    plt.savefig(name + "_tSNE.pdf",dpi=600)
    plt.close("all")

    

def KDE(dataframe,bandwidth = "scott"):
    starttime = datetime.datetime.now()
    data = dataframe
    nbins = np.ceil(data.max())-np.floor(data.min())
    meshgrid = np.meshgrid(*[np.linspace(np.floor(data.iloc[:,i].min()), np.ceil(data.iloc[:,i].max()), int(nbins[i])+1) for i in range(data.shape[1]) ] )
    mgrid = [i.T for i in meshgrid]
    flatten = [i.flatten() for i in mgrid]
    positions = np.vstack(flatten)
    values = data.values.T
    kernel = stats.gaussian_kde(values,bw_method=bandwidth)
    Density = kernel(positions).T
    Density = Density/Density.sum()
    Coordinate = np.array(flatten)
    KernelDensity = pd.DataFrame(np.vstack([Coordinate,Density]).T)
    KernelDensity.columns = [("Dim-" + str(i+1)) for i in range(data.shape[1])] + ["Density"]
    #KernelDensity.to_csv('KernelDensity.txt',sep="\t",index=False,header=True)
    endtime = datetime.datetime.now()
    #print("\nGaussian Kernel Density Estimation Running Time:")
    #print((endtime-starttime))  
    return(KernelDensity)

def sklearnKDE_1(dataframe):
    from sklearn.neighbors import KernelDensity
    starttime = datetime.datetime.now()
    data = dataframe
    nbins = np.ceil(data.max())-np.floor(data.min())
    nbin = nbins[0]+1
    coor = np.linspace(0,10,11)
    b = data.shape[0]
    dim  = data.shape[1]

    COOR = np.array(list(itertools.product(*[coor for i in range(dim)])))
    params = {'bandwidth': np.linspace(0.01, 1, 100)}
    grid = GridSearchCV(KernelDensity(), params)
    grid.fit(data)
    print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))
    kde = grid.best_estimator_
    bandwidth = grid.best_estimator_.bandwidth

    DataRound = deepcopy(dataframe)
    DataRound = np.round(DataRound,0)
    COOR = Neighbors(DataRound.values)

    # kde = KernelDensity(kernel='gaussian', bandwidth=0.01)
    # bandwidth = 0.01

    kde.fit(data)
    log_pdf = kde.score_samples(COOR)
    pdf = 10**log_pdf    
    pdf = pdf/pdf.sum()    
    KernelDensity = pd.DataFrame(COOR)
    KernelDensity["Density"] = pdf
    KernelDensity.columns = [("Dim-" + str(i+1)) for i in range(data.shape[1])] + ["Density"]
    
    cutoff = KernelDensity.Density.max()/data.shape[0]
    KernelDensity =  KernelDensity.loc[KernelDensity.Density >= cutoff,:]
    KernelDensity.index = range(KernelDensity.shape[0])
    # KernelDensity.to_csv('KernelDensity.txt',sep="\t",index=False,header=False)
    endtime = datetime.datetime.now()
    print("\nGaussian Kernel Density Estimation Running Time:")
    print((endtime-starttime))  
    return KernelDensity,COOR,bandwidth
    
def sklearnKDE(dataframe):
    from sklearn.neighbors import KernelDensity
    starttime = datetime.datetime.now()
    data = dataframe
    nbins = np.ceil(data.max())-np.floor(data.min())
    nbin = nbins[0]+1
    coor = np.linspace(0,10,11)
    b = data.shape[0]
    dim  = data.shape[1]

    #COOR = np.array(list(itertools.product(*[coor for i in range(dim)])))
    params = {'bandwidth': np.linspace(0.01, 2, 100)}
    grid = GridSearchCV(KernelDensity(), params)
    grid.fit(data)
    print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))
    kde = grid.best_estimator_
    bandwidth = grid.best_estimator_.bandwidth
    DataRound = deepcopy(dataframe)
    DataRound = np.round(DataRound,0)
    COOR = Neighbors(DataRound.values)
    
    
    #kde = KernelDensity(kernel='gaussian', bandwidth=0.1)
    kde.fit(data)
    log_pdf = kde.score_samples(COOR)
    pdf = 10**log_pdf
    
    pdf = pdf/pdf.sum()
    
    KernelDensity = pd.DataFrame(COOR)
    KernelDensity["Density"] = pdf
    KernelDensity.columns = [("Dim-" + str(i+1)) for i in range(data.shape[1])] + ["Density"]
    
    # cutoff = KernelDensity.Density.max()/data.shape[0]
    # KernelDensity =  KernelDensity.loc[KernelDensity.Density >= cutoff,:]
    # KernelDensity.index = range(KernelDensity.shape[0])
    # KernelDensity.to_csv('KernelDensity.txt',sep="\t",index=False,header=False)
    endtime = datetime.datetime.now()
    print("\nGaussian Kernel Density Estimation Running Time:")
    print((endtime-starttime))  
    return KernelDensity,COOR,bandwidth

def sklearnKDE_fast(dataframe,COOR,bandwidth):
    from sklearn.neighbors import KernelDensity
    starttime = datetime.datetime.now()
    data=dataframe
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth)
    kde.fit(data)
    log_pdf = kde.score_samples(COOR)
    pdf = 10**log_pdf
    
    pdf = pdf/pdf.sum()
    
    KernelDensity = pd.DataFrame(COOR)
    KernelDensity["Density"] = pdf
    KernelDensity.columns = [("Dim-" + str(i+1)) for i in range(data.shape[1])] + ["Density"]
    
    # cutoff = KernelDensity.Density.max()/data.shape[0]
    # KernelDensity =  KernelDensity.loc[KernelDensity.Density >= cutoff,:]
    # KernelDensity.index = range(KernelDensity.shape[0])
    # KernelDensity.to_csv('KernelDensity.txt',sep="\t",index=False,header=False)
    endtime = datetime.datetime.now()
    #print("\nGaussian Kernel Density Estimation Running Time:")
    #print((endtime-starttime))  
    return KernelDensity


def fastKDE(dataframe,nbins = 20):
    starttime = datetime.datetime.now()
    data = dataframe
    kde = FFTKDE(kernel='gaussian', bw = 1)
    x, y = kde.fit(data.values).evaluate(nbins)
    x = pd.DataFrame(x)
    x = x - x.min()
    x = x/x.max() * (nbins -1)
    y = pd.DataFrame(y)
    KernelDensity = pd.concat([x,y],axis=1)
    KernelDensity.columns = [("Dim-" + str(i+1)) for i in range(data.shape[1])] + ["KDE_Density"]
    KernelDensity.to_csv('KernelDensity.txt',sep="\t",index=False,header=False)
    print("\nGaussian Kernel Density Estimation Running Time:")
    print((endtime-starttime))  
    return(KernelDensity)
    

def Rescale(data):
    data = np.floor(data)
    data = data - data.min(axis=0)
    rescaled = deepcopy(data)
    for i in range(data.shape[1]):
        unscaled = data.iloc[:,i].values

        scaled = deepcopy(unscaled)
        index = 0
        for n in range(int(unscaled.max())+1):
            scaled[unscaled == n] = index
            if sum(unscaled == n) > 0:
                index = index + 1
        rescaled.iloc[:,i] = scaled     
    return(rescaled)


def ColourMAP(Cluster):
    Labels = deepcopy(Cluster)
    Uniq = list(set(Labels))
    Uniq.sort()   
    Colors = [plt.cm.Spectral(col) for col in np.linspace(0, 1, len(Uniq))]
    ColorMap = Labels
    for i in range(len(Labels)):
        if Labels[i] == -1:
            ColorMap[i] = (0,0,0,0.75)
        else:
            ColorMap[i] = Colors[Labels[i]]
    del Labels
    return(ColorMap)



def ClusteringEvaluation(PredictLabel,TrueLabel):
    PLabel = pd.read_csv(PredictLabel, header=0, index_col=0, sep=",")
    TLabel = pd.read_csv(TrueLabel, header=0, sep="\s+")
    TLabel.index = ["Cell"+str(i+1) for i in range(1000)]
    ClassifiableRate = PLabel.shape[0]/TLabel.shape[0]
    print("The Classifiable Rate is " + str(ClassifiableRate))
    STLabel = TLabel.loc[PLabel.index,:]
    
    labels_true = STLabel.index
    labels_pred = PLabel.index
    
    name = PredictLabel[0:len(PredictLabel)-4]
    ARI = adjusted_rand_score(labels_true, labels_pred)
    CP  = completeness_score(labels_true, labels_pred)
    FMI = fowlkes_mallows_score(labels_true, labels_pred, sparse=False)
    Homogeneity = homogeneity_score(labels_true, labels_pred)
    NMI = normalized_mutual_info_score(labels_true, labels_pred, average_method='warn')
    AMI = adjusted_mutual_info_score(labels_true, labels_pred, average_method='warn')
    V_Measure = v_measure_score(labels_true, labels_pred, beta=1.0)
    
    Output = open("Clustering_Evaluation.txt",'a')
    line = name +"\t"+ str(ClassifiableRate) +"\t"+ str(ARI) +"\t"+ str(CP) +"\t"+ str(FMI) +"\t"+ str(Homogeneity) +"\t"+ str(AMI) +"\t"+ str(NMI) +"\t"+ str(V_Measure) +"\n"
    Output.write(line)
    Output.close()



def scDensity2D1(data):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    ax.plot(x, y, z, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.5,markersize=4 )
    
    ax.set_xlabel('Dim-1')
    ax.set_ylabel('Dim-2')
    ax.set_zlabel('Density')
    ax.set_title("2D scDensity")
    plt.show()



def scSurface2D1(data):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True) 
    ax.set_xlabel('Dim-1')
    ax.set_ylabel('Dim-2')
    ax.set_zlabel('Density')
    ax.set_title("2D scDensity")
    ax.elev = 45
    ax.azim = 45
    plt.show()

def scDensity2D(data,name,i,Label):
    unique_labels = list(set(Label))
    unique_labels.sort()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for n, col in zip(unique_labels, colors):
        if n == -1:
            col = [0, 0, 0, 1]
        class_mask = [i == n for i in Label]
        data_mask = data.loc[class_mask,:]

        x = data_mask.iloc[:,0]
        y = data_mask.iloc[:,1]
        z = data_mask.iloc[:,2]
        ax.plot(x, y, z, "o", markerfacecolor=tuple(col), markeredgecolor='k',markeredgewidth=0.5,markersize=4 )
     
 
    ax.set_xlabel('Dim-1')
    ax.set_ylabel('Dim-2')
    ax.set_zlabel('Density')
    ax.set_title("2D scDensity")   
    plt.savefig("Density_" + str(i) + ".png",dpi=300)
    plt.close("all")

    

def Density2D(DimensionReduction,name,i):
    # Define custom color map
    startcolor = '#FFFFFF' 
    midcolor = 	'#008B8B'     
    endcolor = '#8B0000'      
    mycmap = color.LinearSegmentedColormap.from_list('MyCMAP',[startcolor,midcolor,endcolor])
    cm.register_cmap(cmap=mycmap)

    # decide the grid number
    x = DimensionReduction.iloc[:,0]
    y = DimensionReduction.iloc[:,1] 
    
    cells = len(x)  
    nbins = int(np.sqrt(cells))

    # gridding the PCA space
    xy = np.vstack([x,y])
    Max = max(x.max(),y.max())
    Min = min(x.min(),y.min())
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    
    xticks = np.round(xi[:,0],5)
    yticks = np.round(yi[0],5)
    xstep = xticks[1] - xticks[0]
    ystep = yticks[1] - yticks[0]
    Xticks = xticks - 0.5*xstep
    Yticks = yticks - 0.5*ystep
    Xticks = np.append(Xticks,Xticks.max()+xstep)
    Yticks = np.append(Yticks,Yticks.max()+ystep)
    
    # count cells in each grid
    zi = np.zeros((nbins,nbins))
    for n in range(DimensionReduction.shape[0]):
        xpos = np.where((DimensionReduction.iloc[n,0]-Xticks <=0) == True)[0][0]-1
        ypos = np.where((DimensionReduction.iloc[n,1]-Yticks <=0) == True)[0][0]-1
        zi[xpos][ypos] = zi[xpos][ypos] +1

    # decide the minimal acceptable cell number in each grid
    cumulative_prob = 0
    for n in range(cells):
        cumulative_prob  = cumulative_prob  + sum(np.random.binomial(n=cells, p=1/(nbins**2), size=100000)==n)/100000
        if cumulative_prob  >=0.99:
            break
    cutoff = n + 1

    heatmap = deepcopy(zi)
    heatmap[heatmap < cutoff] = 0
    numIsland,IslandLabel,mask = numIslands(heatmap)
    IslandLabel = pd.DataFrame(IslandLabel,columns=["Grid_Pos","Label"])
    cumulative_cells = heatmap.sum()/cells
    object_function = cumulative_cells * numIsland
    
    Grid_Pos = []
    Points = deepcopy(DimensionReduction)
    for n in range(Points.shape[0]):
        Grid_X = np.where((Points.iloc[n,0]-Xticks <=0) == True)[0][0]-1
        Grid_Y = np.where((Points.iloc[n,1]-Yticks <=0) == True)[0][0]-1
        Grid_Pos.append([(Grid_X,Grid_Y)])
    Grid_Pos=pd.DataFrame(Grid_Pos,index=DimensionReduction.index,columns=["Grid_Pos"])
    Points = pd.concat([Points, Grid_Pos],axis=1)
    Points = pd.merge(Points,IslandLabel, how='left', on="Grid_Pos")
    Points = Points.fillna(-1)
    Points.index = DimensionReduction.index
    

    cmap = mycmap
    xpos = xi.flatten()
    ypos = yi.flatten()
    zpos = np.zeros_like (xpos)
    dx = xstep
    dy = ystep
    dz = heatmap.flatten()
    
    max_height = np.max(dz)   
    min_height = np.min(dz)
    rgba = [cmap((k-min_height)/max_height) for k in dz] 
    annotate_text_1 = "Cutoff = " + str(cutoff)
    annotate_text_2 = "Cumulative Cells = " + str(round(cumulative_cells,3))
    annotate_text_3 = "Identified Clusters = " + str(numIsland)
    annotate_text_4 = "Object Function = " + str(object_function)
    
    zi_mask = zi * mask
    LabeledPoints = Points.loc[Points["Label"] != -1,:]
    fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz,color=rgba)
    ax.view_init(45, 80)
    plt.title("Peaks")
    plt.xlabel("PC-"+ str(2*i+1))
    plt.ylabel("PC-"+ str(2*i+2))
    plt.gcf().text(0.8, 0.64, annotate_text_1, fontsize=6)
    plt.gcf().text(0.8, 0.61, annotate_text_2, fontsize=6)
    plt.gcf().text(0.8, 0.58, annotate_text_3, fontsize=6)
    plt.savefig(name + "/" + name +  "_3D-Barploot-PC" + str(2*i+1) + "_" + "PC" + str(2*i+2) + ".png",dpi=300)
    plt.close("all")


    fig = plt.figure()
    ax = fig.add_subplot(221,projection='3d')
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz,color=rgba)
    ax.view_init(45, 80)
    ax.set_xlabel("PC-"+ str(2*i+1))
    ax.set_ylabel("PC-"+ str(2*i+2))
    ax.set_title("Peaks")
    
    ax = fig.add_subplot(222)
    ax.imshow(zi_mask.T,cmap='MyCMAP',aspect="auto",vmin=0,vmax=zi.max(),origin="lower")
    #ax.set_xlabel("PC-"+ str(2*i+1))
    ax.set_ylabel("PC-"+ str(2*i+2))
    ax.set_title("Density")

    ax = fig.add_subplot(223)
    ax.plot(x,y,'o', markerfacecolor='darkcyan', markeredgecolor='k', markersize=6) 
    ax.set_xlabel("PC-"+ str(2*i+1))
    ax.set_ylabel("PC-"+ str(2*i+2))
    ax.set_title("All Cells")
    
    ax = fig.add_subplot(224)
    ax.plot(LabeledPoints.iloc[:,0],LabeledPoints.iloc[:,1],'o', markerfacecolor='crimson', markeredgecolor='k', markersize=6)
    ax.set_xlabel("PC-"+ str(2*i+1))
    ax.set_ylabel("PC-"+ str(2*i+2))
    ax.set_title("Filtered Cells")
    ax.set_xlim((min(DimensionReduction.iloc[:,0]),max(DimensionReduction.iloc[:,0])))
    ax.set_ylim((min(DimensionReduction.iloc[:,1]),max(DimensionReduction.iloc[:,1])))

    fig = plt.gcf()
    fig.set_size_inches(10,7)
    plt.suptitle("Clustering Decision Plot on 2 Principle Component")
    plt.savefig(name + "/" + name + "DecisionPlot-PC" + str(2*i+1) + "_" + "PC" + str(2*i+2) + ".png",dpi=300)
    fig.clf()
    plt.close("all")
    
    Decision = np.array(Points["Label"])
    return(numIsland,Decision)



def ElbowPlot(S6,name):
    Dimension = 100
    Dimension = min(Dimension,S6.shape[0])
    plt.figure(figsize=(4,3),dpi=600)
    pca = PCA(n_components=Dimension,)
    PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
    plt.plot(pca.explained_variance_,'o', markerfacecolor='darkcyan',markeredgecolor='k',markeredgewidth=0.5,markersize=3)
    plt.xlabel("Principle Components")
    plt.ylabel("Variance")
    plt.title(name)
    plt.savefig(name+"/ElbowPlot.pdf",dpi=600)
    plt.savefig(name+"/ElbowPlot.tiff",dpi=600)
    plt.close("all")




def CumulativePCAPlot(S6,name):
    Dimension = 50
    pca = PCA(n_components=Dimension)
    PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
    ratio = pca.explained_variance_ratio_
    vector = []
    plt.figure(figsize=(4,3),dpi=600)
    for i in range(len(ratio)):
        vector.append(sum(ratio[0:i]))
    plt.plot(vector,'o', markerfacecolor='darkcyan',markeredgecolor='k',markeredgewidth=0.5,markersize=3)
    plt.xlabel("Principle Components")
    plt.ylabel("Cummulative Explained Variance Ratio")
    plt.title("Variance of each PC")
    plt.savefig(name + "/CummlativePlot.pdf",dpi=600)
    plt.close("all")

    


def PCAHeatmap(S6,name):
    Dimension = 50
    pca = PCA(n_components=Dimension)
    PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
    for i in range(PCA_s.shape[1]):
        Recorder = []
        for j in range(S6.shape[1]):
            corr = np.round(np.corrcoef(PCA_s.iloc[:,i],S6.iloc[:,j])[0][1],4)
            Recorder.append([S6.columns[j],j,corr])
        Recorder = pd.DataFrame(Recorder,columns = ["Gene","Index","Pearson"])
        Recorder = Recorder.sort_values(by="Pearson" , ascending=False)
        h = Recorder.head(25)
        t = Recorder.tail(25)
        recorder = pd.concat([h,t])
        
        Heatmap = S6.iloc[np.random.randint(low=0,high=S6.shape[0],size=500) , recorder.Index]
        
        cm = sns.clustermap(Heatmap.T,xticklabels=False,cmap="viridis",figsize=(7,14))
        cm.ax_row_dendrogram.set_visible(False)
        cm.ax_col_dendrogram.set_visible(False)
        cm.cax.set_visible(False)
        ax = cm.ax_heatmap
        ax.set_xlabel("500 Randomly Selected Cells")
        ax.set_title("Heatmap on Genes Best Correlates with PC-" + str(i+1))
        cm.savefig("PCHeatmap" + str(i+1))


def Scatter_Plot(Data,xlab,ylab):
    plt.figure(figsize=(8,6),dpi=1200)
    plt.plot(Data.iloc[:,0],Data.iloc[:,1],'o', markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.5,markersize=3)
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    
def Labeled_Scatter_Plot(Data,label,xlab,ylab):
    unique_labels = list(set(label))
    unique_labels.sort()
    plt.figure(figsize=(8,6),dpi=1200)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for n, col in zip(unique_labels, colors):
        if n == -1:
            col = [0, 0, 0, 1]
        class_mask = [i == n for i in label]
        Data_mask = Data.loc[class_mask,:]
    plt.plot(Data_mask.iloc[:,0],Data_mask.iloc[:,1],'o', markerfacecolor=tuple(col), markeredgecolor='k',markeredgewidth=0.5,markersize=3)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    


def Silhouette_Plot(PCA_s,name,Label):
    Label = pd.DataFrame(Label)
    Label = Label.loc[Label.iloc[:,0]!=-1,]
    PCA_s = PCA_s.iloc[Label.index,: ]
    PCA_s.index = np.linspace(start=0,stop=PCA_s.shape[0]-1,num=PCA_s.shape[0], dtype=int)
    Label.index = np.linspace(start=0,stop=Label.shape[0]-1,num=Label.shape[0], dtype=int)
    n_clusters = len(set(Label.iloc[:,0]))
    

    silhouette_avg = silhouette_score(PCA_s, Label)
    sample_silhouette_values = silhouette_samples(PCA_s, Label)
    
    fig, ax1 = plt.subplots(1)
    fig.set_size_inches(10,8)
    ax1.set_xlim([-0.1, 1])
    ax1.set_ylim([0, len(PCA_s) + (n_clusters + 1) * 10])

    y_lower = 10
    for i in range(n_clusters):
        ith_cluster_silhouette_values = sample_silhouette_values[ Label.loc[Label.iloc[:,0]==i,].index]
        ith_cluster_silhouette_values.sort()
    
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),0, ith_cluster_silhouette_values, facecolor=color, edgecolor=color, alpha=0.7)
        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
    
    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")
    
    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="brown", linestyle="--")
    ax1.set_yticks([]) 
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    
    
     
    fig = plt.gcf()
    fig.set_size_inches(6,10)
    plt.title("Silhouette Scores")
    plt.savefig("Silhouette.png",dpi=300)
    fig.clf()
    plt.close("all")
    
    

def Shannon_Index(prob):
    prob = prob[prob>0]
    temp = [-i*np.log2(i) for i in prob]
    SI = np.round(sum(temp),4)
    return(SI)


def K_Estimation_Old(PCAData,Kstart=10):
    K = Kstart
    min_eigen = 0
    while min_eigen < 0.0001 :
        #print(K)
        KNNetwork = kneighbors_graph(PCAData, n_neighbors=K, include_self=False)
        KNNeighbor =[set(KNNetwork[i].nonzero()[1]) for i in range(len(PCAData))]
        Ajacency = np.array([[len(KNNeighbor[i].intersection(KNNeighbor[j]))/K for j in range(len(KNNeighbor))] for i in range(len(KNNeighbor))])
        np.fill_diagonal(Ajacency,0)      
        prune = 5/(K+K-5)
        Ajacency[Ajacency < prune] = 0
        Ajacency[Ajacency > 0] = 1
        Degree = np.zeros_like(Ajacency)
        np.fill_diagonal(Degree, Ajacency.sum(axis=1))
        L = Degree - Ajacency
        L_csr = sparse.csr_matrix(L)
        L_eigs = sparse.linalg.eigsh(L_csr, k = L.shape[0]-1, return_eigenvectors=False)
        min_eigen = L_eigs.min()
        #print(min_eigen)
        if min_eigen >= 0.0001:
            break
        else:
            K = K + 1
    print("The estimated K is " + str(K))
    print("Eigenvalue: "+str(min_eigen))
    return(K)
    gc.collect()


def K_Estimation(PCAData,Kstart=10):
    K = Kstart
    while True :
        KNNetwork = kneighbors_graph(PCAData, n_neighbors=K, include_self=False)
        KNNeighbor =[set(KNNetwork[i].nonzero()[1]) for i in range(len(PCAData))]
        starttime = datetime.datetime.now()
        SNN= np.array([[len(KNNeighbor[i].intersection(KNNeighbor[j]))/K for j in range(len(KNNeighbor))] for i in range(len(KNNeighbor))])
        SNN = pd.DataFrame(SNN)    
        SNN = np.triu(SNN,1)
        SNN = pd.DataFrame(SNN)
        Edges = SNN.stack()
        Edges.index.rename(['NodeA', 'NodeB'], inplace=True)
        Edges = Edges.to_frame('Weight').reset_index()
        Edges = Edges.loc[Edges.NodeA != Edges.NodeB,:]   
        prune = 2/(K+K-2)
        Edges = Edges.loc[Edges.Weight >= prune,:]
        Edges.index = np.arange(Edges.shape[0])    
        edge_list = [(Edges.iloc[i,0],Edges.iloc[i,1]) for i in range(Edges.shape[0])]
        G = nx.from_edgelist(edge_list)
        connect = nx.is_connected(G)
        #print(connect)
        if connect == True:
            break
        else:
            K = K + 1
    print("The estimated K is " + str(K))
    return(K)
    gc.collect()

def KNN_Density_Ball_SNN(PCAData,K):
    starttime = datetime.datetime.now()
    N = PCAData.shape[0]
    d = PCAData.shape[1]
    KNN_Dist = kneighbors_graph(PCAData, n_neighbors=K, mode='distance', metric='minkowski', p=2, include_self=False)
    Distance = np.round(squareform(pdist(np.array(PCAData), metric='euclidean')),4)
    R = np.array([])
    for i in range(Distance.shape[0]):
        line = Distance[i]
        line.sort()
        R = np.append(R,line[K])
        
    v = (pi**(d/2))/gamma(d/2+1)*(R**d)
    del Distance,R
    Den = (K/float(N))/v
    Den = Den/Den.sum()
    Den = pd.Series(Den)
    KNNetwork = kneighbors_graph(PCAData, n_neighbors=K, include_self=False)
    KNNeighbor =[set(KNNetwork[i].nonzero()[1]) for i in range(len(PCAData))]
    starttime = datetime.datetime.now()
    SNN= np.array([[len(KNNeighbor[i].intersection(KNNeighbor[j]))/K for j in range(len(KNNeighbor))] for i in range(len(KNNeighbor))])
    SNN = pd.DataFrame(SNN)    
    SNN = np.triu(SNN,1)
    SNN = pd.DataFrame(SNN)
    Edges = SNN.stack()
    Edges.index.rename(['NodeA', 'NodeB'], inplace=True)
    Edges = Edges.to_frame('Weight').reset_index()
    Edges = Edges.loc[Edges.NodeA != Edges.NodeB,:]   
    prune = 2/(K+K-2)
    Edges = Edges.loc[Edges.Weight >= prune,:]
    Edges["NodeADen"] = Den[Edges.NodeA].values
    Edges["NodeBDen"] = Den[Edges.NodeB].values
    ## To define the higher vertice as the "From" vertice, the lower vertice as the "To" vertice 
    Diff = (Edges["NodeADen"] - Edges["NodeBDen"]) < 0
    Edges.loc[Diff,['NodeA','NodeB']] = Edges.loc[Diff,['NodeB','NodeA']].values
    Edges.loc[Diff,['NodeADen','NodeBDen']] = Edges.loc[Diff,['NodeBDen','NodeADen']].values
    Edges["EdgeDen"] = (Edges["NodeADen"] + Edges["NodeBDen"])/2
    Edges = Edges.sort_values(by=["EdgeDen","NodeADen"] , ascending = False)
    Edges.index = np.arange(Edges.shape[0])    
    gc.collect()
    endtime = datetime.datetime.now()
    # print("\nShared Nearest Networt Construction Time:")
    # print((endtime-starttime))  
    return(Edges,Den)


def SNN_Graph(PCA_Data,K):
    starttime = datetime.datetime.now()
    KNNetwork = kneighbors_graph(PCA_Data, n_neighbors=K, include_self=False)
    KNNeighbor =[set(KNNetwork[i].nonzero()[1]) for i in range(len(PCA_Data))]
    starttime = datetime.datetime.now()
    SNN= np.array([[len(KNNeighbor[i].intersection(KNNeighbor[j]))/K for j in range(len(KNNeighbor))] for i in range(len(KNNeighbor))])
    SNN = pd.DataFrame(SNN)
        

    Den = (SNN.sum(axis=1)-1)/2
    Den = Den/Den.sum()*100000
    SNN = np.triu(SNN,1)
    
    SNN = pd.DataFrame(SNN)
    Edges = SNN.stack()
    Edges.index.rename(['NodeA', 'NodeB'], inplace=True)
    Edges = Edges.to_frame('Weight').reset_index()
    Edges = Edges.loc[Edges.NodeA != Edges.NodeB,:]   
    prune = 2/(K+K-2)
    #prune = 0.6
    Edges = Edges.loc[Edges.Weight >= prune,:]
    Edges["NodeADen"] = Den[Edges.NodeA].values
    Edges["NodeBDen"] = Den[Edges.NodeB].values
    ## To define the higher vertice as the "From" vertice, the lower vertice as the "To" vertice 
    Diff = (Edges["NodeADen"] - Edges["NodeBDen"]) < 0
    Edges.loc[Diff,['NodeA','NodeB']] = Edges.loc[Diff,['NodeB','NodeA']].values
    Edges.loc[Diff,['NodeADen','NodeBDen']] = Edges.loc[Diff,['NodeBDen','NodeADen']].values
    Edges = Edges.sort_values(by="NodeADen" , ascending = False)
    Edges.index = np.arange(Edges.shape[0])
    
    lowest = Edges.NodeBDen.argmin()
    lowestNodeB = int(Edges.values[lowest,1])
    lowestNodeBDen = Edges.values[lowest,4]
    nodes = np.unique(Edges.NodeA)
    np.random.shuffle(nodes)
    nodes = nodes[0:50]
    append = pd.DataFrame(columns=Edges.columns)
    append.NodeA = nodes
    append.NodeB = lowestNodeB
    append.Weight = prune
    append.NodeADen = Den[append.NodeA].values
    append.NodeBDen = lowestNodeBDen
    
    Edges = pd.concat([Edges,append],axis=0,ignore_index=True)
    Edges = Edges.sort_values(by="NodeADen" , ascending = False)
    Edges.index = np.arange(Edges.shape[0])    
    endtime = datetime.datetime.now()
    print("\nShared Nearest Networt Construction Time:")
    print((endtime-starttime))  
    return(Edges,Den)


def superlevel_set_filtration_plot(PCA_Data,Edges,Den,N):
    def root(i):
        while i != idx[i]:
            idx[i] = idx[idx[i]]
            i = idx[i]
        return i
    
    
    def union(f,t):
        F = root(f)
        T = root(t)
        idx[T] = F
    
    colors = ["crimson","coral","dodgerblue","purple","cyan"]
    ID = np.array(list(range(PCA_Data.shape[0])))
    idx = deepcopy(ID)
    #SNN, Den = KNN_Density_Ball_SNN(Data, K = 30)
    data = TSNE(n_components=2,random_state=1).fit_transform(PCA_Data)
    Data = pd.DataFrame(data, columns=["tSNE-1","tSNE-2"]) 
    
    edge_list = [(Edges.iloc[i,0],Edges.iloc[i,1]) for i in range(Edges.shape[0])]
    # Generate coordiante of nodes and edges
    # 3D coordinates of nodes
    node_xyz = np.concatenate([Data.values, Den.values.reshape(-1,1)],axis=1)
    # 3D coordinates of edges
    edge_xyz = np.array([(node_xyz[u], node_xyz[v]) for u, v in edge_list])
    # 2D coordinates of nodes
    node_xy = deepcopy(data)
    # 2D coordinates of edges
    edge_xy = np.array([(node_xy[u], node_xy[v]) for u, v in edge_list])    

    fig = plt.figure(figsize=(12,5),clear=True)
    ## the 3D network
    # the first n edges of highest density
    ax1 = fig.add_subplot(1,3,1, projection="3d")
    edge_highlight = edge_xyz[0:N+1]
    raw_shape = edge_highlight.shape
    npoints = raw_shape[0] * raw_shape[1]
    node_highlight = edge_highlight.reshape(1,npoints,3)[0]
    # Union Find
    edges = edge_list[0:N]
    for i in range(N):
        union(edges[i][0],edges[i][1])
    idx = [root(n) for n in idx ]
    seq, count = np.unique(idx, return_counts=True)
    codes = (np.where(count>1))[0]
    L = len(codes)
    Codes = {}
    for l in range(L):
        Codes[l] = seq[codes[l]]      

    # Plot the nodes - alpha is scaled by "depth" automatically
    ax1.plot(*node_xyz.T, "o", markersize = 6, mec = "black", mew = 0.3, mfc="lightgrey")
    # Labeling highlight nodes
    for l in range(L):
        ax1.plot(node_xyz[ID[idx == Codes[l]].tolist(),0], node_xyz[ID[idx == Codes[l]].tolist(),1], node_xyz[ID[idx == Codes[l]].tolist(),2], "o", markersize = 6, mec = "black", mew = 0.3, mfc=colors[l])
    xticks =  np.linspace(data.min(axis=0)[0], data.max(axis=0)[0], 10)
    yticks =  np.linspace(data.min(axis=0)[1], data.max(axis=0)[1], 10)
    xx, yy = np.meshgrid(xticks, yticks)
    den = node_highlight.min(axis=0)[2]
    zz = np.ones_like(xx) * den
    ax1.plot_surface(xx, yy, zz, alpha=0.2, color = "springgreen")
    # Formatting axis
    ax1.zaxis.set_rotate_label(False) 
    ax1.set_xlabel("Dimension-1", fontsize=16, labelpad=10)
    ax1.set_ylabel("Dimension-2", fontsize=16, labelpad=10)
    ax1.set_zlabel("Density",fontsize=16, labelpad=10, rotation=90)
    ax1.view_init(30,30)
    
    ## 2D network
    ax2 = fig.add_subplot(1,3,2)
    # Plot the nodes - alpha is scaled by "depth" automatically
    edge_highlight = edge_xy[0:N+1]
    raw_shape = edge_highlight.shape
    npoints = raw_shape[0] * raw_shape[1]
    node_highlight = edge_highlight.reshape(1,npoints,2)[0]
    ax2.plot(data[:,1], data[:,0], "o", markersize = 6, mec = "black", mew = 0.3, mfc="lightgrey")
    # Labeling highlight nodes
    for l in range(L):
        ax2.plot(node_xy[ID[idx == Codes[l]].tolist(),1], node_xy[ID[idx == Codes[l]].tolist(),0], "o", markersize = 6, mec = "black", mew = 0.3, mfc=colors[l])
    # Formatting axis
    ax2.set_xlabel("Dimension-2", fontsize=16, labelpad=10)
    ax2.set_ylabel("Dimension-1", fontsize=16, labelpad=10)
    ax2.set_xlim(data.min(axis=0)[1],data.max(axis=0)[1])
    ax2.set_ylim(data.min(axis=0)[0],data.max(axis=0)[0])
    ax2.invert_yaxis() 

    ## 2D Connected component
    ax3 = fig.add_subplot(1,3,3)
    edge_highlight = edge_xy[0:N+1]
    # Labeling highlight nodes
    for l in range(L):
        ax3.plot(node_xy[ID[idx == Codes[l]].tolist(),1], node_xy[ID[idx == Codes[l]].tolist(),0], "o", markersize = 6, mec = "black", mew = 0.3, mfc=colors[l])
    # Plot the edges
    for edge in edge_highlight:
        ax3.plot(edge[:,1], edge[:,0], linestyle="dashed", color="lightgrey", linewidth=0.2)
    # Formatting axis
    ax3.set_xlabel("Dimension-2", fontsize=16, labelpad=10)
    ax3.set_ylabel("Dimension-1", fontsize=16, labelpad=10)
    ax3.set_xlim(data.min(axis=0)[1],data.max(axis=0)[1])
    ax3.set_ylim(data.min(axis=0)[0],data.max(axis=0)[0])
    ax3.invert_yaxis() 

    # ax4 = fig.add_subplot(1,4,4)
    # Bar_Len = []
    # for l in range(L):
    #     bar_id = ID[idx == Codes[l]].tolist()
    #     bar_max = Den[bar_id].max()
    #     bar_min = Den[bar_id].min()
    #     bar_len = bar_max - bar_min
    #     Bar_Len.append(bar_len)
    
    # ax4.barh(list(np.arange(L)), Bar_Len, 0.2, color="w")
    # #ax4.barh(["1","2","3"], [Den[n1].max()-Den[n1].min(), Den[n2].max()-Den[n2].min(), Den[n3].max()-Den[n3].min()], 0.2, left=[Den[n1].min(), Den[n2].min(), Den[n3].min()], color=["crimson","coral","dodgerblue"])
    # ax4.invert_yaxis() 
    # ax4.invert_xaxis() 
    # ax4.ticklabel_format(style='scientific', scilimits=(0,0), useMathText=True, axis="x")
    # ax4.set_xlabel("Barcode Length", fontsize=16, labelpad=10)
    # ax4.set_ylabel("Persistence Barcodes ($H_{0}$)", fontsize=16, labelpad=10)
    # ax4.set_xlim(Den.max(),Den.min())


    
    ### Overall settings
    fig.tight_layout()#调整整体空白
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)  #调整子图间距
    plt.savefig("Superlevel_Set_Filtration_" + str(N) + "_edges.png", dpi=600)
    plt.close("all")



def PersistentHomology(Data, tSNE, Edges, Den, name, pd_name ="raw", Stratified = True, iter_plot = False, cluster_plot = True):
    def root(i):
        while i != ID[i]:
            ID[i] = ID[ID[i]]
            i = ID[i]
        return i

    def connected(ID,f,t):
        return(root(ID,f) == root(ID,t))
    
    def union(f,t):
        F = root(f)
        T = root(t)
        ID[T] = F

    def deunion(x):
        ID[x] = x

    starttime = datetime.datetime.now()
    #tsne = TSNE(n_components=2, random_state=1).fit_transform(Data)
    #tsne_embedded = pd.DataFrame(tsne,columns=["tSNE-1","tSNE-2"]) 
    
    ID = np.arange(Data.shape[0])
    size = np.ones_like(ID,dtype="int")
    visit = -1*np.ones_like(ID,dtype="int")
    
    BirthIndex = []
    DeathIndex = []
    BirthDensity = []
    DeathDensity = []
    PeakSize = []
    CurrentLevelCluster = []
    UpperLeverCluster = []
    Barcode=[]
    barcode = 0
    Cluster = {0:np.ones_like(ID,dtype="int")*Den.values.argmax()}
    PeakSize.append(len(Den))
    Cluster_U = deepcopy(Cluster)
    Label = np.zeros_like(ID)
    CurrentLevelCluster.append(Den.argmax())
    UpperLeverCluster.append("None")
    BirthIndex.append(Den.argmax())
    DeathIndex.append(Den.argmin())
    BirthDensity.append(Den.max())
    DeathDensity.append(Den.min())
    Barcode.append(barcode)
    first = 0

    
    for i in range(Edges.shape[0]):    
        f = Edges.iloc[i,0]
        t = Edges.iloc[i,1]
        F = root(f)
        T = root(t)
        visit[f] = 1
        visit[t] = 1
        
        if Den[F] > Den[T]:
            union(F,T)
            size[F] += size[T]
            if Den[T] != Den[t]:
                deunion(T)      
                if F == Den.argmax() and first == 0:
                    first = first + 1
                    Label_0 = np.array([root(n) for n in ID ])
                    Label_0[Label_0 != F] = 0
                    Cluster_U[0] = Label_0
                    
                Label = np.array([root(n) for n in ID ])
                Label[visit < 0] = -1    
                if iter_plot == True:
                    superlevel_set_filtration_plot(Data,Edges,Den,i)
                BirthIndex.append(T)
                DeathIndex.append(t)
                BirthDensity.append(Den[T])
                DeathDensity.append(Den[t])
                PeakSize.append(size[T])
                barcode += 1
                Barcode.append(barcode)
                Label_C = deepcopy(Label)
                Label_U = deepcopy(Label)
                Label_C[Label_C != T] = 0
                Label_U[Label_U != F] = 0
                Cluster[barcode] = Label_C
                Cluster_U[barcode] = Label_U
                UpperLeverCluster.append(F)
                CurrentLevelCluster.append(T)
                union(F,T)


        if Den[F] < Den[T]:
            union(T,F)
            size[T] += size[F]
            if Den[T] != Den[t]:
                deunion(F)
                if T == Den.argmax() and first == 0:
                    first = first + 1
                    Label_0 = np.array([root(n) for n in ID ])
                    Label_0[Label_0 != T] = 0
                    Cluster_U[0] = Label_0        
                
                Label = np.array([root(n) for n in ID ])
                Label[visit < 0] = -1   
                if iter_plot == True:
                    superlevel_set_filtration_plot(Data,Edges,Den,i)
                BirthIndex.append(F)
                DeathIndex.append(t)
                BirthDensity.append(Den[F])
                DeathDensity.append(Den[t])  
                PeakSize.append(size[F])
                barcode += 1
                Barcode.append(barcode)     
                Label_C = deepcopy(Label)
                Label_U = deepcopy(Label)
                Label_C[Label_C != F] = 0
                Label_U[Label_U != T] = 0
                Cluster[barcode] = Label_C
                Cluster_U[barcode] = Label_U
                UpperLeverCluster.append(T)
                CurrentLevelCluster.append(F)
                union(T,F)
                
        if Den[F] == Den[T]:
            union(F,T)
    #print(Label)
    
    if barcode >= 1:
        FinalLabel = np.array([root(n) for n in ID ])
        unknownPeakIndex = list(set(FinalLabel) - set(Label))
        knownPeakIndex = list(set(Label))
        if -1 in knownPeakIndex:
            knownPeakIndex.remove(-1)
        if len(unknownPeakIndex) >= 1:
            for ind in unknownPeakIndex:
                barcode = barcode + 1
                Barcode.append(barcode)
                birth_ind = Den[FinalLabel == ind].argmax()
                death_ind = Den[FinalLabel == ind].argmin()
                birth_den = Den[FinalLabel == ind].max()
                death_den = Den[FinalLabel == ind].min()
                cluster_size = sum(FinalLabel == ind)
                
                Label_C = deepcopy(FinalLabel)
                Label_C[Label_C != ind] = 0
                Label_U = deepcopy(FinalLabel)
                BirthIndex.append(birth_ind)
                DeathIndex.append(death_ind)
                BirthDensity.append(birth_den)
                DeathDensity.append(death_den)
                PeakSize.append(cluster_size)     
                Cluster[barcode] = Label_C
                Cluster_U[barcode] = Label_U
                UpperLeverCluster.append(T)
                CurrentLevelCluster.append(F)        


    PH = pd.DataFrame({"Barcode":Barcode,"BirthInd":BirthIndex,"DeathInd":DeathIndex,"BirthDen":BirthDensity,"DeathDen":DeathDensity,"CurrentCluster":CurrentLevelCluster,"UpperCluster":UpperLeverCluster,"PeakSize":PeakSize})
    PH["BarcodeLength"] = PH["BirthDen"] - PH["DeathDen"]
    PH = PH.sort_values(by="PeakSize" , ascending=False)
    # cluster0id = PH.loc[0,"CurrentCluster"]
    # upperind = PH.loc[(PH.UpperCluster == PH.loc[0,"CurrentCluster"]),"Barcode"].values[0]
    # Cluster0 = Cluster_U[upperind]
    # Cluster0[Cluster0 != cluster0id ]= 0
    # Cluster[0] = Cluster0
    # Cluster = pd.DataFrame(Cluster)
    


    if Stratified == True:
        nClust = PH.shape[0]
        Cluster = {key: value for key, value in Cluster.items() if key in PH.index}
        font1 = {'color':  'crimson', 'weight': 'normal'}
        font2 = {'color':  'lightgrey', 'weight': 'normal'}
        idx = PH["PeakSize"].idxmax()
        Maximum =  PH.BirthDen.max()
        step = Maximum/100
        plt.figure(figsize=(6,4),dpi=600)
        plt.plot([0, Maximum], [0, Maximum], ls="--", color="black")
        for i in list(PH.index)[::-1]:
            if PH.PeakSize[i] >= 10:
                plt.plot(PH.DeathDen,PH.BirthDen, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.1, markersize=3 )
                plt.text(PH.DeathDen[i]+step, PH.BirthDen[i]+step, r'${:}$'.format(PH.PeakSize[i]), fontdict=font1)
            else:
                plt.plot(PH.DeathDen,PH.BirthDen, "o", markerfacecolor='lightgrey', markeredgecolor='k',markeredgewidth=0.1, markersize=3 )
                plt.text(PH.DeathDen[i]+step, PH.BirthDen[i]+step, r'${:}$'.format(PH.PeakSize[i]), fontdict=font2)
        plt.title("Persistent Diagram, PeaksNum = " + str(nClust))
        plt.ylabel("Birth")
        plt.xlabel("Death")
        plt.xticks()
        plt.yticks()
        plt.savefig(name +"/Persistent_Diagram_" + str(pd_name) + "_stratified.pdf", dpi=600)
        #plt.savefig(name +"/Persistent_Diagram_" + str(pd_name) + "_stratified.tiff", dpi=600)
        plt.close("all")
        
        
        # PH = PH.iloc[1:PH.shape[0],:]
        # Cluster = {key: value for key, value in Cluster.items() if key in PH.index}
        # font1 = {'color':  'crimson', 'weight': 'normal'}
        # font2 = {'color':  'lightgrey', 'weight': 'normal'}
        # idx = PH["PeakSize"].idxmax()
        # Maximum =  PH.BirthDen.max()
        # step = Maximum/100
        # plt.figure(figsize=(6,4),dpi=600)
        # plt.plot([0, Maximum], [0, Maximum], ls="--", color="black")
        # for i in list(PH.index)[::-1]:
        #     if PH.PeakSize[i] >= 10:
        #         plt.plot(PH.DeathDen,PH.BirthDen, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.1, markersize=3 )
        #         plt.text(PH.DeathDen[i]+step, PH.BirthDen[i]+step, r'${:}$'.format(PH.PeakSize[i]), fontdict=font1)
        #     else:
        #         plt.plot(PH.DeathDen,PH.BirthDen, "o", markerfacecolor='lightgrey', markeredgecolor='k',markeredgewidth=0.1, markersize=3 )
        #         plt.text(PH.DeathDen[i]+step, PH.BirthDen[i]+step, r'${:}$'.format(PH.PeakSize[i]), fontdict=font2)
        # plt.title("Persistent Diagram, PeaksNum = " + str(nClust))
        # plt.ylabel("Birth")
        # plt.xlabel("Death")
        # plt.xticks()
        # plt.yticks()
        # plt.savefig(name +"/Persistent_Diagram_" + str(pd_name) + "_stratified-2.pdf", dpi=600)
        # #plt.savefig(name +"/Persistent_Diagram_" + str(pd_name) + "_stratified.tiff", dpi=600)
        # plt.close("all")

    else:

        font1 = {'color':  'red', 'weight': 'normal'}
        font2 = {'color':  'lightgrey', 'weight': 'normal'}
        idx = PH["PeakSize"].idxmax()
        Maximum =  PH.BirthDen.max()
        step = Maximum/100
        plt.figure(figsize=(6,4),dpi=600)
        plt.plot([0, Maximum], [0, Maximum], ls="--", color="black")
        plt.plot(PH.DeathDen,PH.BirthDen, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.1, markersize=3 )
        for i in list(PH.index)[::-1]:
            plt.text(PH.DeathDen[i]+step, PH.BirthDen[i]+step, r'${:}$'.format(PH.PeakSize[i]), fontdict=font1)
        plt.title("Persistent Diagram, PeaksNum = " + str(PH.shape[0]))
        plt.ylabel("Birth")
        plt.xlabel("Death")
        plt.xticks()
        plt.yticks()
        plt.savefig(name +"/Persistent_Diagram_" + str(pd_name) + ".pdf", dpi=600)
        plt.savefig(name +"/Persistent_Diagram_" + str(pd_name) + ".tiff", dpi=600)
        plt.close("all")


    if cluster_plot == True:
        colors = ["firebrick","darkcyan","purple","darkslategrey","salmon","palegreen","darkorange","lightseagreen","blue","coral","chocolate","tan","darkgoldenrod","olivedrab","yellow","seagreen","steelblue","royalblue","slateblue","violet","crimson", "firebrick","darkcyan","purple","darkslategrey","salmon","palegreen","darkorange","lightseagreen","blue","coral","chocolate","tan","darkgoldenrod","olivedrab","yellow","seagreen","steelblue","royalblue","slateblue","violet","crimson"]
        for key in list(Cluster.keys()):
            selected = (Cluster[key]!=0)
            if sum(selected) > 1:
                plt.figure(figsize=(6,4),dpi=600)
                plt.plot(tSNE.iloc[:,1],tSNE.iloc[:,0],'o', markerfacecolor="lightgrey", markeredgecolor='k',markeredgewidth=0.1,markersize=3)
                plt.plot(tSNE.iloc[selected,1],tSNE.iloc[selected,0],'o', markerfacecolor=colors[key], markeredgecolor='k',markeredgewidth=0.1,markersize=3)
                plt.title("Cluster " + str(key))
                plt.xlabel("tSNE-2")
                plt.ylabel("tSNE-1")
                ax = plt.gca()
                ax.invert_yaxis() 
                plt.savefig(name + "/tSNE_runid=" + str(pd_name) + "_clusterid=" + str(key) + ".pdf", dpi=600)
                plt.savefig(name + "/tSNE_runid=" + str(pd_name) + "_clusterid=" + str(key) + ".tiff", dpi=600)
                plt.close("all")

    
    endtime = datetime.datetime.now()
    # print("\nTopological Data Analysis Running Time:")
    # print((endtime-starttime)) 
    return PH, Cluster, Cluster_U



def DownSamplingStability(PCA_Data, PH, Cluster, CellID, tSNE, K, name, iteration = 10):
      
    ## KNN_Density_Ball_SNN
    ## probability landscape assumption
    ID = np.arange(PCA_Data.shape[0])
    frac = 0.9
    Index_sampling = {}
    PD_sampling = {}
    Cluster_sampling = {}
    nCells = PCA_Data.shape[0]
    print("Downsampling Times: ")
    for run in np.arange(iteration):
        print(run)
        sample_ind = list(np.random.choice(list(CellID.index), size=int(nCells*frac), replace=False) )
        #sample_ind = list(np.random.choice(list(CellID.index), size=int(nCells*frac), replace=False) )
        sample_ind.sort()
        subsample  = PCA_Data.iloc[sample_ind,:]
        subCellID  = CellID.iloc[sample_ind,:]
        subCellID.index = np.arange(int(nCells*frac))  
        subtSNE = tSNE.iloc[sample_ind,:]
        subtSNE.index = subCellID.index
        subEdges,subDen = KNN_Density_Ball_SNN(subsample, K)
        #subEdges, subDen = SNN_Graph(subsample, K)
        subDen.index = subCellID.index
        subPH,subCluster,subCluster_U = PersistentHomology(subsample, subtSNE, subEdges, subDen, name, "sampling_iter="+str(run), Stratified = False, iter_plot = False, cluster_plot = False)
        subPH = subPH.sort_values(by="PeakSize" , ascending=False)
        Index_sampling[run] = sample_ind
        Selected_subBarcode = subPH.loc[subPH.PeakSize >= 5,"Barcode"].values
        subCluster = {key: value for key, value in subCluster.items() if key in Selected_subBarcode}
        
        PD_sampling[run] = subPH
        Cluster_sampling[run] = subCluster
        gc.collect()
        #print("****************************************************************************************************************")
 
    
    table=[]
    Selected_Barcode = PH.loc[PH.PeakSize >= 10,"Barcode"].values
    selectedPH = PH.loc[PH.PeakSize >= 10,:]
    Cluster = {key: value for key, value in Cluster.items() if key in Selected_Barcode}
    for ref in Cluster.keys():
        jaccardindex = []
        RefCluster = Cluster[ref]
        for run in Cluster_sampling.keys():
            index = []
            for query in Cluster_sampling[run].keys():
                Refcluster = RefCluster[Index_sampling[run]]
                Querycluster = Cluster_sampling[run][query]
                intersection = sum((Refcluster*Querycluster) > 0)
                union = sum((Refcluster+Querycluster) > 0)
                if union == 0:
                    jaccard = 0
                else:
                    jaccard = np.round(intersection/union,4)
                index.append(jaccard)
            index = np.array(index)
            Jaccard = index.max()
            query = index.argmax()
            jaccardindex.append(Jaccard)
            table.append((ref,run,query,Jaccard))           
    table = pd.DataFrame(table)
    table.columns = ["Ref","Run","Subsample","Jaccard"]
    table.to_csv(name + "/Jaccard.csv")
    
    summary = pd.DataFrame(table.groupby('Ref').describe()["Jaccard"])
    summary["cv"] = summary["std"]/summary["mean"]
    summary = pd.merge(selectedPH,summary,left_index=True,right_index=True,how='outer')
    summary.to_csv(name + "/Jaccard_Distribution.csv")
    
    plt.figure(figsize=(12, 3), dpi=600)
    sns.boxplot(x="Ref", y="Jaccard", data=table)
    #JaccardTable.boxplot()
    # for i in range(selectedPH .shape[0]):
    #     y = table.loc[table.Ref == selectedPH.iloc[i,0], ["Jaccard"]]
    #     x = np.random.normal(i, 0.06, size=len(y))
    #     plt.plot(x, y, 'o', markerfacecolor='crimson',
    #              markeredgecolor='black', markeredgewidth=0.3, markersize=6)
    plt.xlabel("Cluster")
    plt.ylabel("Stability")
    plt.xticks()
    plt.yticks()
    plt.ylim(0, 1)
    plt.savefig(name + "/Downsampling_stability.pdf", dpi=600)
    plt.close("all")

    res = summary.loc[:, ["Barcode", "CurrentCluster",
                          "PeakSize", "mean", "std", "min", "50%", "max", "cv"]]
    res["Stable_Fraction"] = res["PeakSize"] * res["mean"]
    res = res.sort_values(by="Stable_Fraction" , ascending=False)
    #res.to_csv(name + "/Jaccard_Distribution.csv")
    # fig = plt.figure(figsize=(18,6),dpi=600)
    # plt.bar(list(summary.index), summary["mean"], yerr=summary["std"], align='center', alpha=0.5, ecolor='black', capsize=10)
    # plt.ylabel('Jaccard Index')
    # plt.ylim(0,1)
    # plt.xlabel("Cluster", fontsize=16, labelpad=10)
    # plt.ylabel("Stability", fontsize=16, labelpad=10)
    # plt.xticks(list(summary.index),fontsize = 16)
    # plt.xticks(fontsize=16)
    # plt.yticks(fontsize=16)
    # # Save the figure and show
    # plt.tight_layout()
    # plt.savefig(name + "_downsampling_stability_errbar.png", dpi=600)
    # plt.show()
    print(res)
    return res





def Find_Daughter(pylist,Tree):
    L = np.unique(pylist)
    if len(L) == 1:
            return 0
    else:
        dellist = np.ones_like(pylist)
        for i in range(pylist.shape[0]):
            try:
                if Tree[pylist[i]] in pylist:
                    dellist[i] = 0
            finally:
                continue
        return(pylist[dellist>0][0])
        



def Peaks_Cells_Label(PH,Cluster,Cluster_U,candidatePeaks):
    U = list(set(PH.UpperCluster))
    U.remove("None")
    
    ## For peaks other than peaks0
    U0 = PH.CurrentCluster[0]
    if len(U) >= 1:
        U.remove(U0)
    for i in range(len(U)):
        #print(i)
        parent_ind = PH.loc[PH.CurrentCluster==U[i],["Barcode"]]
        parent_ind = parent_ind.values.flatten()
        if parent_ind.shape[0] > 0:
            parent_ind = parent_ind[0]
            daughter_ind = PH.loc[PH.UpperCluster==U[i],["Barcode"]]
            daughter_ind = daughter_ind.values.flatten()
            # if len(daughter_ind) > 1:
            #     print(i)
            daughter_ind_min = np.array([sum(Cluster_U[i]>0) for i in daughter_ind]).argmin()
            daughter_ind_min = daughter_ind[daughter_ind_min]
            Cluster[parent_ind] = Cluster_U[daughter_ind_min]
            # daughter_ind_max = np.array([sum(Cluster_U[i]>0) for i in daughter_ind]).argmax()
            # daughter_ind_max = daughter_ind[daughter_ind_max]
            # Cluster[parent_ind] = Cluster_U[daughter_ind_max]
    ## For peak0
    Cluster[0] = Cluster_U[0]
    
    ## Summary 2 results
    Cluster = pd.DataFrame(Cluster)
    PH_Peaks = PH.loc[candidatePeaks,:]
    Cluster_Peaks = Cluster.loc[:,PH_Peaks.index]
    
    Tree = dict(zip(list(PH_Peaks.UpperCluster)+[0], list(PH_Peaks.CurrentCluster)+[0]))
    if Cluster_Peaks.shape[1] == 1:
        Cluster_Peaks["Label"] = Cluster_Peaks.iloc[:,0]
    else:
        Cluster_Peaks["Label"] =  [Find_Daughter(Cluster_Peaks.values[i],Tree) for i in range(Cluster_Peaks.shape[0])] 
    return PH_Peaks,Cluster_Peaks


def Rename_Peaks(Cluster_Peaks):
    uniq_label = list(np.unique(Cluster_Peaks.Label.values))

    if len(uniq_label) == 1 and uniq_label[0] !=0:
        rename_dict = dict(zip(uniq_label, [1]))
    else:
        rename_dict = dict(zip(uniq_label, np.arange(len(uniq_label))))
    Cluster_Peaks.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks.Label.values]
    
    return Cluster_Peaks    
    
def Iterative_KNN_Calling0(cluster_peaks, PCA_Data, tSNE, name):
    Cluster_Peaks = deepcopy(cluster_peaks)
    Prediction = Cluster_Peaks.Label.values
    if sum(Prediction==0)==0:
        pass
    else:
        Prediction[Prediction==0] = -1
        tSNEPlot(tSNE, name +"/Label_of_Peaks", Prediction)
        All_List = list(Cluster_Peaks.index)
        UnKnown_List = list(np.where(Prediction<0)[0])
        Known_List  = list(np.where(Prediction>0)[0])
        while len(UnKnown_List) >0 :
            UnKnown = PCA_Data.iloc[UnKnown_List,:]
            Known   = PCA_Data.iloc[Known_List,:]
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(Known)
            dist, index = nbrs.kneighbors(UnKnown)
            dist  = dist.flatten()
            index = index.flatten()
            known_nearest_index = index[dist.argmin()]
            new_nearest_index = UnKnown_List[dist.argmin()]
            Prediction[new_nearest_index] = Prediction[Known_List[known_nearest_index]] 
            UnKnown_List = list(np.where(Prediction<0)[0])
            Known_List  = list(np.where(Prediction>0)[0])
            if len(UnKnown_List)%100 == 0:
                print(len(UnKnown_List))
                tSNEPlot(tSNE, name + "/tSNE_Unknown=" + str(len(UnKnown_List)),Prediction) 
        Cluster_Peaks["Label"] = Prediction
    tSNEPlot(tSNE, name + "/Label_of_AllCells",Cluster_Peaks["Label"])
    return Cluster_Peaks


def Iterative_KMeans_Calling(cluster_peaks, PCA_Data, tSNE, name):
    Cluster_Peaks = deepcopy(cluster_peaks)
    Prediction = Cluster_Peaks.Label.values
    if sum(Prediction==0)==0:
        pass
    else:
        Prediction[Prediction==0] = -1
        tSNEPlot(tSNE, name +"/Label_of_Peaks", Prediction)
        All_List = np.array(Cluster_Peaks.index)
        UnKnown_List = np.where(Prediction<0)[0]
        Known_List  =  np.where(Prediction>0)[0]
        while len(UnKnown_List) >0 :
            UnKnown = PCA_Data.iloc[UnKnown_List,:]
            Known   = PCA_Data.iloc[Known_List,:]
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(Known)
            dist, index = nbrs.kneighbors(UnKnown)
            index = index[:,0]
            Prediction[UnKnown_List] = Prediction[Known_List[index]] 
            UnKnown_List = list(np.where(Prediction<0)[0])
            Known_List  = list(np.where(Prediction>0)[0])
            if len(UnKnown_List)%100 == 0:
                print(len(UnKnown_List))
                tSNEPlot(tSNE, name + "/tSNE_Unknown=" + str(len(UnKnown_List)),Prediction) 
        Cluster_Peaks["Label"] = Prediction
    tSNEPlot(tSNE, name + "/Label_of_AllCells",Cluster_Peaks["Label"])
    return Cluster_Peaks

def Iterative_KNN_Calling(cluster_peaks, PCA_Data, tSNE, name, iter_plot = False):
    Cluster_Peaks = deepcopy(cluster_peaks)
    Prediction = Cluster_Peaks.Label.values
    if sum(Prediction==0)==0:
        pass
    else:
        Prediction[Prediction==0] = -1
        tSNEPlot(tSNE, name +"/Label_of_Peaks", Prediction)
        All_List = np.array(Cluster_Peaks.index)
        UnKnown_List = np.where(Prediction<0)[0]
        Known_List  =  np.where(Prediction>0)[0]
        while len(UnKnown_List) >0 :
            UnKnown = PCA_Data.iloc[UnKnown_List,:]
            Known   = PCA_Data.iloc[Known_List,:]
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(Known)
            dist, index = nbrs.kneighbors(UnKnown)
            dist  = dist.flatten()
            index = index.flatten()
            step_len = 100
            ind = np.argsort(dist)[0:step_len]
            known_nearest_index = index[ind]
            new_nearest_index = UnKnown_List[ind]
            Prediction[new_nearest_index] = Prediction[Known_List[known_nearest_index]] 
            UnKnown_List = np.where(Prediction<0)[0]
            Known_List  =  np.where(Prediction>0)[0]
            print(len(UnKnown_List))
            if iter_plot == True:
                tSNEPlot(tSNE, name + "/tSNE_Unknown=" + str(len(UnKnown_List)),Prediction) 
        Cluster_Peaks["Label"] = Prediction
    tSNEPlot(tSNE, name + "/Label_of_AllCells",Cluster_Peaks["Label"])
    return Cluster_Peaks




def scRNA_Seq_TDA(file):
    starttime = datetime.datetime.now()
    ### Step-1 Specify dataset
    name = file.split("_")[0] 
    if not os.path.exists(name):
        os.makedirs(name)
         
    ### Step-2 Data processing
    S1 = pd.read_csv(file, header=0, index_col=0, sep=",")
    gene = name  + "_gene.csv"
    Gene = pd.read_csv(gene, header=0, index_col=False, sep=",")
    Symbol = dict(zip(Gene.iloc[:,0],Gene.iloc[:,1]))
    S1.index = [Symbol[i] for i in S1.index]
    #S2 = Quality_Control(S1,cellNum=0,UMINumLowLimit=500,UMINumUpLimit=10000,mtPct=0.2)
    S3 = Normalize_Rawdata(S1)
    S4 = Select_HEGenes(S3,ngenes=2500)
    S5 = Scaling_Data(S4)
    S6 = S5.T   
    TrueLabel = name  + "_barcodes_anno.csv"
    TrueLabel = pd.read_csv(TrueLabel, header=0, index_col=False, sep=",")
    truelabel = TrueLabel.values[:,1]
    
    
    ### Step-3 PCA dimension reduction for density peaks finding
    Dimension = 100
    pca = PCA(n_components=Dimension,  random_state = 1) 
    PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
    ElbowPlot(S6,name)
    Dim = PCA_Permutation(S6, name, nperm = 100, nPC = 40)
    
    # The estimated number of PCs is 12
    PCA_Data = np.round(PCA_s.iloc[:,0:Dim],4)
    CellID = pd.DataFrame(PCA_Data.index, columns=["ID"])
    nCells = PCA_Data.shape[0]
    
    
    ### Step-4 tSNE dimension reduction for visualization
    tsne = TSNE(n_components=2, random_state=1).fit_transform(PCA_Data)
    tSNE = pd.DataFrame(tsne,columns=["tSNE-1","tSNE-2"]) 
    tSNEPlot(tSNE, name + "/TrueLabel_of_AllCells",truelabel,legend={"loc":'lower left', "bbox_to_anchor":(1,-0.03),"fontsize":5, "ncol":1})
    endtime = datetime.datetime.now()
    print("\nPreprocessing Time:")
    print((endtime-starttime))  
    
    
    ### Step-5 Node-weighted SNN Graph
    starttime = datetime.datetime.now()
    K = K_Estimation(PCA_Data,10)
    endtime = datetime.datetime.now()
    print("\nNeighborhood Size Esitmation Time:")
    print((endtime-starttime))  
    
    # the least K that makes SNN a 1-connected graph is estimated at 20
    starttime = datetime.datetime.now()
    Edges, Den = KNN_Density_Ball_SNN(PCA_Data, K)
    #Edges, Den = SNN_Graph(PCA_Data, K)
    endtime = datetime.datetime.now()
    print("\nSNN network construction Time:")
    print((endtime-starttime))  
    
    
    ### Step-6 Persistenct homology on superlevel set filtration
    starttime = datetime.datetime.now()
    PH, Cluster, Cluster_U = PersistentHomology(PCA_Data, tSNE, Edges, Den, name, pd_name = "raw", Stratified = False, iter_plot = False, cluster_plot = False)   
    endtime = datetime.datetime.now()
    print("\nPeaks Identification Time:")
    print((endtime-starttime))  
    
    ### Step-7 
    # Optional: Cluster stability evaluation via downsampling 
    starttime = datetime.datetime.now()
    Summary = DownSamplingStability(PCA_Data, PH, Cluster, CellID, tSNE, K, name, iteration = 20)
    N = sum(Summary.Estimation>=10)   
    endtime = datetime.datetime.now()
    print("\nPeaks Identification Time:")
    print((endtime-starttime))  
    # OR: Specify the Number of Clusters N from PeakSize of PH table
    # N = (defined by user, learnt from PH table) 
    
    
    ### Step-8 Report peaks_fraction of each cluster
    starttime = datetime.datetime.now()
    PH_Peaks, Cluster_Peaks = Peaks_Cells_Label(PH,Cluster,Cluster_U,selectPeaks)
    uniq_label = list(np.unique(Cluster_Peaks.Label.values))

    ### Step-9 Rename peaks
    rename_dict = dict(zip(uniq_label, np.arange(len(uniq_label))))
    Cluster_Peaks.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks.Label.values]
    Cluster_Peaks.to_csv(name + "/ClusteringLabelPeaks.csv")

    ### Step-9 Iterative KNN to call all cells back
    Cluster_All = Iterative_KNN_Calling(Cluster_Peaks, PCA_Data, tSNE, name, iter_plot=False)
    endtime = datetime.datetime.now()
    print("\nAssign All Other Cells Time:")
    print((endtime-starttime)) 
    
    ### Step-10 Output results
    #ARI = np.round(adjusted_rand_score(truelabel,Cluster_ALL.Label),4)
    starttime = datetime.datetime.now()
    g=open(name + "/summary.txt","w")
    g.write("The estimated number of PCs is: " + str(Dim))
    g.write("The estimated neighborhood is: " + str(K))
    g.write("The estimated cluster number is: " + str(N))
    #g.write(str(ARI))
    g.close()

    PCA_Data["RealLabel"] = truelabel
    PCA_Data["PeakLabel"] = Cluster_Peaks.Label.values
    PCA_Data["AllLabel"] = Cluster_All.Label.values
    PCA_Data.to_csv(name + "/PCA.csv")
    
    S4 = S4.T
    S4["RealLabel"] = truelabel
    S4["PeakLabel"] = Cluster_Peaks.Label.values
    S4["AllLabel"]  = Cluster_All.Label.values
    S4.to_csv(name + "/S4.csv")
    
    tSNE.to_csv(name + "/tSNE.csv")
    PH.to_csv(name + "/PersistenceDiagram.csv")
    PH_Peaks.to_csv(name + "/PersistenceDiagramPeaks.csv")
    Cluster_All.to_csv(name + "/ClusteringLabel.csv")
    PredictLabel = Cluster_All.Label.values
    LabelID = np.unique(PredictLabel)
    RawData = S1.T


    ### Step-11 Output subcluster in raw UMI countmatrix format
    for id in LabelID:
        print(id)
        subRawData = RawData.iloc[PredictLabel == id,:]
        subRawData.T.to_csv(name + "/subCluster_" + str(id) + ".txt", sep="\t")
        subtSNE = tSNE.iloc[PredictLabel == id,:]
        subtSNE["Label"] = 0
        subtSNE.to_csv(name + "/tSNE_" + str(id) + ".txt", sep="\t")
    print("The estimated number of PCs is: " + str(Dim))
    print("The estimated neighborhood is: " + str(K))
    print("The estimated cluster number is: " + str(N))
    endtime = datetime.datetime.now()
    print("\nOutput Results Time:")
    print((endtime-starttime)) 
    return Cluster_All,PH

def scRNA_Seq_TDA_sub(subid):
    
    ### Step-1 Specify dataset
    file = "subCluster_" + str(subid) + ".txt"
    tsne = "tSNE_" + str(subid) + ".txt"
    name = str(subid)
    if not os.path.exists(name):
        os.makedirs(name)
    
    ### Step-2 Data processing
    S1 = pd.read_csv(file, header=0, index_col=0, sep="\s+")
    #S2 = Quality_Control(S1)
    S3 = Normalize_Rawdata(S1)
    S4 = Select_HEGenes(S3,ngenes=2500)
    S5 = Scaling_Data(S4)
    S6 = S5.T   
    
    
    ### Step-3 PCA dimension reduction for density peaks finding
    Dimension = 100
    Dimension = min(Dimension,S6.shape[0])
    pca = PCA(n_components=Dimension,  random_state = 1) 
    PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
    ElbowPlot(S6,name)
    Dim = PCA_Permutation(S6, name, nperm = 100)

    
    # The estimated number of PCs is 12
    PCA_Data = np.round(PCA_s.iloc[:,0:Dim],4)
    CellID = pd.DataFrame(PCA_Data.index, columns=["ID"])
    nCells = PCA_Data.shape[0]
    
    
    ### Step-4 tSNE dimension reduction for visualization
    ### Step-4 tSNE dimension reduction for visualization
    tSNE = pd.read_csv(tsne, header=0, index_col=0, sep="\s+")
    tSNEPlot(tSNE, name + "/subcluster_" +str(subid),tSNE.Label)
    
    
    ### Step-5 Node-weighted SNN Graph
    K = K_Estimation(PCA_Data,5)

    # the least K that makes SNN a 1-connected graph is estimated at 20
    Edges, Den = KNN_Density_Ball_SNN(PCA_Data, K)
    #Edges, Den = SNN_Graph(PCA_Data, K)
    
    ### Step-6 Persistenct homology on superlevel set filtration
    PH, Cluster, Cluster_U = PersistentHomology(PCA_Data, tSNE, Edges, Den, name, pd_name = "raw", Stratified = False, iter_plot = False, cluster_plot = False)
    
    
    ### Step-7 
    # Optional: Cluster stability evaluation via downsampling 
    Summary = DownSamplingStability(PCA_Data, PH, Cluster, CellID, tSNE, K, name, iteration = 20)
    N = sum(Summary.Estimation>=10)   
    # OR: Specify the Number of Clusters N from PeakSize of PH table
    # N = (defined by user, learnt from PH table) 
    
    
    ### Step-8 Report peaks_fraction of each cluster
    PH_Peaks, Cluster_Peaks = Peaks_Cells_Label(PH,Cluster,Cluster_U,selectPeaks)
    uniq_label = list(np.unique(Cluster_Peaks.Label.values))

    ### Step-9 Rename peaks
    rename_dict = dict(zip(uniq_label, np.arange(len(uniq_label))))
    Cluster_Peaks.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks.Label.values]
    Cluster_Peaks.to_csv(name + "/ClusteringLabelPeaks.csv")
    

    ### Step-10 Iterative KNN to call all cells back
    Cluster_All = Iterative_KNN_Calling(Cluster_Peaks, PCA_Data, tSNE, name)

    ### Step-11 Output results
    #ARI = np.round(adjusted_rand_score(truelabel,Cluster_ALL.Label),4)
    g=open(name + "/summary.txt","w")
    g.write("The estimated neighborhood is: " + str(K) +"\n")
    g.write("The estimated number of PCs is: " + str(Dim) +"\n")
    g.write("The estimated cluster number is: " + str(N) +"\n")
    #g.write(str(ARI))
    g.close()
    
    PCA_Data["PeakLabel"] = Cluster_Peaks.Label.values
    PCA_Data["AllLabel"] = Cluster_All.Label.values
    PCA_Data.to_csv(name + "/PCA.csv")
    tSNE.to_csv(name + "/tSNE.csv")
    PH.to_csv(name + "/PersistenceDiagram.csv")
    PH_Peaks.to_csv(name + "/PersistenceDiagramPeaks.csv")
    Cluster_All.to_csv(name + "/ClusteringLabel.csv")
    PredictLabel = Cluster_All.Label.values
    LabelID = np.unique(PredictLabel)
    RawData = S1.T

    ### Step-11 Output subcluster in raw UMI countmatrix format
    for id in LabelID:
        print(id)
        subRawData = RawData.iloc[PredictLabel == id,:]
        subRawData.T.to_csv(name + "/subCluster_" + str(id) + ".txt", sep="\t")
        subtSNE = tSNE.iloc[PredictLabel == id,:]
        subtSNE["Label"] = 0
        subtSNE.to_csv(name + "/tSNE_" + str(id) + ".txt", sep="\t")
    print("The estimated number of PCs is: " + str(Dim))
    print("The estimated neighborhood is: " + str(K))
    print("The estimated cluster number is: " + str(N))


def IntraClusterCorrelation(PCA_Data,Label_Column,name):
    """
    Parameters
    ----------
    PCA_Data : pandas dataframe
        PCA dimension reduced data with group labels
    Label_Column : str
        "RealLabel": ground-truth cell type label
        "PeakLabel": TDA-predicted cell type labels of peak cells
        "AllLabel":  TDA-predicted cell type labels of all cells 

    Returns
    -------
    None.

    """
    
    ID = list(set(PCA_Data.loc[:, Label_Column].values.tolist()))
    if 0 in ID:
        ID.remove(0)

    CORR = pd.DataFrame(columns=["Corr", "Label"])
    datacolumns = [c for c in PCA_Data.columns if c.startswith(
        "PC")] + [Label_Column]
    Data = PCA_Data.loc[:, datacolumns]
    for i in ID:
        print(i)
        cluster = Data.loc[Data.loc[:, Label_Column] == i, :]
        cluster = cluster.iloc[:, 0:cluster.shape[1]-1]
        corr = np.corrcoef(cluster)
        if cluster.shape[0] > 1 :
            corr = np.triu(corr, 1)
        corr = corr.reshape(1, -1)
        corr = corr[corr > 0]
        corr = pd.DataFrame(corr, columns=["Corr"])
        corr["Label"] = i
        CORR = pd.concat([CORR, corr], axis=0, sort=False)

    plt.figure(figsize=(12, 3), dpi=600)
    sns.boxplot(x="Label", y="Corr", data=CORR)
    plt.xlabel("Cluster")
    plt.ylabel("Intra-cluster pairwise correlation")
    plt.xticks(rotation=90)
    plt.yticks()
    plt.ylim(0, 1)
    plt.savefig(name+"/Intracluster_Consistency_" + Label_Column + ".pdf", dpi=600)
    plt.close("all")


def Heatmap(PCA_Data, Label_Column, legend={"loc":'center left', "bbox_to_anchor":(-0.1,0.5),"fontsize":7, "ncol":1}):
    data = deepcopy(PCA_Data)
    data = data.loc[data.loc[:,Label_Column]!=0,:]
    columns = list(data.columns)
    data = data.sort_values(by=Label_Column , ascending=True)
    species = data.loc[:,Label_Column]
    unique_labels = list(set(species))
    unique_labels.sort()
    Colors = [plt.cm.Spectral(col) for col in np.linspace(0, 1, len(unique_labels))]
    lut = dict(zip(species.unique(), Colors))
    row_colors = species.map(lut)
    data = data.loc[:,[i for i in columns if i.startswith("PC")]]
    
    plt.figure(figsize=(4,4),dpi=600)
    g = sns.clustermap(data, col_cluster = False, row_cluster = False, yticklabels=False, row_colors=row_colors,edgecolor="white",  cmap='seismic',  z_score=1, cbar_pos=(1.05, 0.1, 0.03, 0.65))
    for label in unique_labels:
        g.ax_row_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(loc=legend["loc"], bbox_to_anchor=legend["bbox_to_anchor"],fontsize=legend["fontsize"], ncol=legend["ncol"])
    
    plt.savefig("Heatmap_" + Label_Column + ".pdf", dpi=600)
    plt.close("all")



def ClusterSpecificGene(S4,threshold = 0.5):
    data = deepcopy(S4)
    data = data.loc[data.loc[:,"PeakLabel"]!=0,:]
    columns = list(data.columns)
    data = data.sort_values(by="PeakLabel" , ascending=True)
    mat = data.groupby(['PeakLabel']).sum()
    mat = np.round(mat/mat.sum(axis=0),4)
    mat = mat.T
    gene = np.array(mat.index)
    marker = dict()
    for i in range(mat.shape[1]):
        values = mat.iloc[:,i]
        values = values.sort_values(ascending=False)
        marker[i+1] = np.array(values[values>=threshold].index)
    return mat, marker


def ClusterSpecificGeneSet(S4, name, top=10, group=None, pct_cutoff = 0.1, log2FC_cutoff = 0.25, pvalue_adjust_cutoff = 0.05 ):
    data = deepcopy(S4)
    data = data.sort_values(by="PeakLabel" , ascending=True)
    gene = [i for i in data.columns.to_list() if not i.endswith("Label") ]
    Gene = data.loc[:,gene]
    counts = data.PeakLabel.value_counts()
    
    # if some cluster has less than 10 cells, assign neareast unknown cells to it 
    # for better differential expression genes calculation
    small  = counts[counts<10]
    if len(small) > 0:
        for c in small.index.to_list():
            known = Gene.loc[data.PeakLabel == c,:]
            unknown = Gene.loc[data.PeakLabel == 0,:]
            dist = cdist(known,unknown).min(axis=0)
            unknown_index= np.argsort(dist)
            data_index   = unknown.iloc[unknown_index,:].index.to_list()
            temp = data.loc[data_index,:]
            temp["PeakLabel"] = c
            temp = temp.loc[temp.PeakLabel == temp.AllLabel,:]
            temp = temp.head(10-small[c])
            add_index = temp.index.to_list()
            data.loc[add_index,["PeakLabel"]] = c
    
    # selects peak cells for differential genes finding
    data = data.loc[data.loc[:,"PeakLabel"]!=0,:]
    data = data.sort_values(by="PeakLabel", ascending=True)
    cluster_label = np.unique(data.PeakLabel.values)
    
    DiffGene = {}
    DataFrame = pd.DataFrame()
    
    if group == None:
        for clust in cluster_label:
            print("Finding differentially expressed genes for cluster " + str(clust))
            result = pd.DataFrame(None,columns=["Cluster","Gene","log2FC","lgo2FC.abs","PCT.1","PCT.2","pvalue","pvalue_adjust"])
            log2FC = []
            PCT1  = []
            PCT2  = []
            pvalue = []
            pvalue_adjust = []
            temp = deepcopy(data)
            cluster = temp.loc[temp.PeakLabel == clust,:]
            others  = temp.loc[temp.PeakLabel != clust,:]
            for g in gene:
                c = cluster[g].values
                o = others[g].values
                #stat = stats.ranksums(c,o,alternative="two-sided")
                stat = stats.ranksums(c,o,alternative="greater")
                log2FC.append(np.round(np.log2((c.mean()+0.0001)/(o.mean()+0.001)),3))
                PCT1.append(np.round(sum(c>0)/len(c),3))
                PCT2.append(np.round(sum(o>0)/len(o),3))
                pvalue.append(stat.pvalue)
            pvalue_adjust = multipletests(pvalue, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[1]
            result["Cluster"] = [clust]*len(gene)
            result["Gene"] = gene
            result["log2FC"] = log2FC
            result["log2FC.abs"] = np.abs(log2FC)
            result["PCT.1"] = PCT1 
            result["PCT.2"] = PCT2
            result["pvalue"] = pvalue
            result["pvalue_adjust"] = pvalue_adjust
            result = result.sort_values(by="log2FC.abs" , ascending=False)
            result = result.loc[np.abs(result["log2FC"]) >= log2FC_cutoff, :]
            result = result.loc[(result["PCT.1"]  >= pct_cutoff) + (result["PCT.2"] >= pct_cutoff), :]
            #result = result.loc[np.abs(result["PCT.2"])  >= pct_cutoff , :]
            result = result.loc[result["pvalue_adjust"]  <= pvalue_adjust_cutoff , :]
            result.index = np.arange(result.shape[0])
            DiffGene[clust] = result
            DataFrame = pd.concat([DataFrame,result],axis=0)
        DataFrame.index = np.arange(DataFrame.shape[0])
        DataFrame.to_csv(name+"/DifferentialGeneSets.csv")
        
        # select top 10 diff genes for each cluster
        High_Diff_Gene = []
        for clust in DiffGene:
            High_Diff_Gene.extend(list(DiffGene[clust]["Gene"].values[0:top]))
            #High_Diff_Gene.extend(list(DiffGene[clust]["Gene"]))
            
        # average transcriptome for each cluster
        Gene = data.loc[:, gene]
        Node = pd.DataFrame(index=gene)
        for clust in cluster_label:
            select = (data.PeakLabel == clust)
            select_ind = select[select == True].index.to_list()
            cluster = Gene.loc[select_ind, :]
            cluster_mean = pd.DataFrame(cluster.mean(axis=0), columns=[clust])
            Node = pd.concat([Node, cluster_mean], axis=1)
        Node.to_csv(name+"/Cluster_Transcriptome.csv")
        # plot dendrogram

        #Dendrogram = Node.loc[High_Diff_Gene,:]
        Gene = Gene.T
        Dendrogram = Gene.loc[High_Diff_Gene,:]
        # plt.figure(figsize=(8,6),dpi=300)
        # fig, ax = plt.subplots(1, 1)
        # g = sns.clustermap(Gene, center = 4,cmap="viridis",
        #                    col_cluster = True, row_cluster = False,
        #                    col_colors = col_colors,
        #                    dendrogram_ratio=(.1, .2),
        #                    cbar_pos=(.02, .32, .03, .2),figsize=(6,6))
        # #g.ax_col_dendrogram.remove()
        # g.ax_heatmap.set_xticks([])
        # Col_Colors = np.array(col_colors)
        # borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(Col_Colors)])
        # for b0, b1, label in zip(borders[:-1], borders[1:], unique_labels):
        #     g.ax_col_colors.text((b0 + b1) / 2, 1.06, label, color='black', ha='center', va='bottom', transform=g.ax_col_colors.get_xaxis_transform())
        # plt.savefig("Hierarchical_Dendrogram_clustering.pdf", dpi=600)
        # plt.close("all")
        
        # species = pd.Series(Node.columns.to_list())
        # unique_labels = list(set(species))
        # unique_labels.sort()
        # Colors = [plt.cm.Spectral(col) for col in np.linspace(0, 1, len(unique_labels))]
        # lut = dict(zip(species, Colors))
        # col_colors = species.map(lut)
        # plt.figure(figsize=(8,6),dpi=300)
        # fig, ax = plt.subplots(1, 1)
        # g = sns.clustermap(Node, center = 4,cmap="viridis",
        #                    col_cluster = True, row_cluster = False,
        #                    col_colors = col_colors,
        #                    dendrogram_ratio=(.1, .2),
        #                    cbar_pos=(.02, .32, .03, .2),figsize=(6,6))
        # #g.ax_col_dendrogram.remove()
        # g.ax_heatmap.set_xticks([])
        # Col_Colors = np.array(col_colors)
        # borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(Col_Colors)])
        # for b0, b1, label in zip(borders[:-1], borders[1:], unique_labels):
        #     g.ax_col_colors.text((b0 + b1) / 2, 1.06, label, color='black', ha='center', va='bottom', transform=g.ax_col_colors.get_xaxis_transform())
        # plt.savefig("Hierarchical_Dendrogram_clustering-2.pdf", dpi=600)
        # plt.close("all")
    
        species = data.PeakLabel
        unique_labels = list(set(species))
        unique_labels.sort()
        Colors = [plt.cm.Spectral(col) for col in np.linspace(0, 1, len(unique_labels))]
        lut = dict(zip(species.unique(), Colors))
        col_colors = species.map(lut)
        
        plt.figure(figsize=(12,64),dpi=300)
        fig, ax = plt.subplots(1, 1)
        g = sns.clustermap(Dendrogram, center = 4,cmap="viridis",
                           col_cluster = False, row_cluster = False,
                           col_colors = col_colors,
                           dendrogram_ratio=(.1, .2),
                           cbar_pos=(.02, .32, .03, .2),figsize=(12,64))
        g.ax_row_dendrogram.remove()
        g.ax_col_dendrogram.remove()
        g.ax_heatmap.set_xticks([])
        Col_Colors = np.array(col_colors)
        borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(Col_Colors)])
        for b0, b1, label in zip(borders[:-1], borders[1:], unique_labels):
            g.ax_col_colors.text((b0 + b1) / 2, 1.06, label, color='black', ha='center', va='bottom', transform=g.ax_col_colors.get_xaxis_transform())
        plt.savefig(name + "/Hierarchical_Dendrogram_sequential.pdf", dpi=600)
        plt.close("all")
        
                
        plt.figure(figsize=(8,6),dpi=300)
        fig, ax = plt.subplots(1, 1)
        g = sns.clustermap(Dendrogram, center = 4,cmap="viridis",
                           col_cluster = False, row_cluster = False,
                           col_colors = col_colors,
                           dendrogram_ratio=(.1, .2),
                           cbar_pos=(.02, .32, .03, .2),figsize=(8,6))
        g.ax_row_dendrogram.remove()
        g.ax_col_dendrogram.remove()
        g.ax_heatmap.set_xticks([])
        Col_Colors = np.array(col_colors)
        borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(Col_Colors)])
        for b0, b1, label in zip(borders[:-1], borders[1:], unique_labels):
            g.ax_col_colors.text((b0 + b1) / 2, 1.06, label, color='black', ha='center', va='bottom', transform=g.ax_col_colors.get_xaxis_transform())
        plt.savefig(name+"/Hierarchical_Dendrogram_sequentia2.pdf", dpi=600)
        plt.close("all")
    
    else:
        C = group[0]
        O = group[1]
        print("Finding differentially expressed genes between cluster " + str(C) + " and cluster " + str(O) )
        result = pd.DataFrame(None,columns=["Gene","log2FC","log2FC.abs","PCT.1","PCT.2","pvalue","pvalue_adjust"])
        log2FC = []
        PCT1  = []
        PCT2  = []
        pvalue = []
        pvalue_adjust = []
        temp = deepcopy(data)
        cluster_c = temp.loc[temp.PeakLabel == C,:]
        cluster_o = temp.loc[temp.PeakLabel == O,:]
        for g in gene:
            c = cluster_c[g].values
            o = cluster_o[g].values
            stat = stats.ranksums(c,o,alternative="two-sided")
            log2FC.append(np.round(np.log2((c.mean()+0.0001)/(o.mean()+0.001)),3))
            PCT1.append(np.round(sum(c>0)/len(c),3))
            PCT2.append(np.round(sum(o>0)/len(o),3))
            pvalue.append(stat.pvalue)
        pvalue_adjust = multipletests(pvalue, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[1]
        result["Gene"] = gene
        result["log2FC"] = log2FC
        result["log2FC.abs"] = np.abs(log2FC)
        result["PCT.1"] = PCT1 
        result["PCT.2"] = PCT2
        result["pvalue"] = pvalue
        result["pvalue_adjust"] = pvalue_adjust
        result = result.sort_values(by="log2FC.abs" , ascending=False)
        # result = result.loc[np.abs(result["log2FC"]) >= log2FC_cutoff, :]
        # result = result.loc[np.abs(result["PCT.1"])  >= pct_cutoff , :]
        # result = result.loc[np.abs(result["PCT.2"])  >= pct_cutoff , :]
        # result = result.loc[np.abs(result["pvalue_adjust"])  <= pvalue_adjust_cutoff , :]
        result.index = np.arange(result.shape[0])
        result.to_csv(name+"/DifferentialGeneBetweenCluster_" + str(C) +"-"+ str(O) +".csv")
    return DiffGene


def tSNEforMarker0(tSNE, S4, gene):
    data = deepcopy(tSNE)
    data["Color"] = S4.loc[:,gene].values
    data = data.sort_values("Color",ascending=True)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgrey","crimson"])
    plt.figure(figsize=(4,3),dpi=600)
    plt.scatter(data.iloc[:,0],data.iloc[:,1], marker='o', s=2, c=data["Color"], cmap = cmap )
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    plt.title(gene)
    plt.savefig("tSNE_" + gene + ".pdf",dpi=600)
    plt.close("all")

def tSNEforMarker(tSNE, S4, gene, cluster):
    data = deepcopy(tSNE)
    data["Color"] = S4.loc[:,gene].values
    data = data.sort_values("Color",ascending=True)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgrey","crimson"])
    plt.figure(figsize=(4,3),dpi=600)
    plt.scatter(data.iloc[:,0],data.iloc[:,1], marker='o', s=2, c=data["Color"], cmap = cmap )
    plt.xlabel("tSNE-1")
    plt.ylabel("tSNE-2")
    plt.title("Cluster_" + str(cluster) + ": " + gene)
    plt.savefig("tSNE_Cluster=" + str(cluster) + "_"  + gene + ".pdf",dpi=600)
    plt.close("all")
