#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 17:24:28 2022

@author: mac
"""

### Step-0 Import packages
import os
os.chdir("/Users/mac/Desktop/scRNA/RealData/MCA/Adult-Stomach/")
from script.functions import *


### Step-1 Specify dataset
starttime = datetime.datetime.now()
file = "Adult-Stomach_dge.csv"
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
# The estimated number of PCs is 30
Dim = 30
PCA_Data = np.round(PCA_s.iloc[:, 0:Dim], 4)
CellID = pd.DataFrame(PCA_Data.index, columns=["ID"])
nCells = PCA_Data.shape[0]


### Step-4 tSNE dimension reduction for visualization
tsne = TSNE(n_components=2, random_state=1).fit_transform(PCA_Data)
tSNE = pd.DataFrame(tsne, columns=["tSNE-1", "tSNE-2"])
tSNEPlot(tSNE, name + "/TrueLabel_of_AllCells", truelabel,
         legend={"loc": 'lower left', "bbox_to_anchor": (1, -0.03), "fontsize": 5, "ncol": 1})
endtime = datetime.datetime.now()
print("\nPreprocessing Time:")
print((endtime-starttime))


### Step-5 Node-weighted SNN Graph
starttime = datetime.datetime.now()
K = K_Estimation(PCA_Data,5)
# The estimated neighborhood size is 10
K = 5
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
PH, Cluster, Cluster_U = PersistentHomology(PCA_Data, tSNE, Edges, Den, name, pd_name = "raw", Stratified = True, iter_plot = False, cluster_plot = False)   
endtime = datetime.datetime.now()
print("\nPeaks Identification Time:")
print((endtime-starttime))  

### Step-7 Optional: Cluster stability evaluation via downsampling 
starttime = datetime.datetime.now()
Summary = DownSamplingStability(PCA_Data, PH, Cluster, CellID, tSNE, K, name, iteration = 20)
candidatePeaks = Summary.loc[Summary.Stable_Fraction >= 10,"Barcode"].values


### Step-8 Report peaks_fraction of each cluster
starttime = datetime.datetime.now()
PH_Peaks, Cluster_Peaks = Peaks_Cells_Label(PH,Cluster,Cluster_U,candidatePeaks)
uniq_label = list(np.unique(Cluster_Peaks.Label.values))
endtime = datetime.datetime.now()
print("\nPeaks Identification Time:")
print((endtime-starttime))  

### Step-9 Rename peaks
rename_dict = dict(zip(uniq_label, np.arange(len(uniq_label))))
Cluster_Peaks.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks.Label.values]


### Step-9 Iterative KNN to call all cells back
Cluster_All = Iterative_KNN_Calling(Cluster_Peaks, PCA_Data, tSNE, name, iter_plot=False)
nClust = len(list(np.unique(Cluster_Peaks.Label.values)))
endtime = datetime.datetime.now()
print("\nAssign All Other Cells Time:")
print((endtime-starttime)) 

### Step-10 Output results
starttime = datetime.datetime.now()
# report parameters
g=open(name + "/summary.txt","w")
g.write("The estimated number of PCs is: " + str(Dim) + "\n")
g.write("The estimated neighborhood is: " + str(K)+ "\n")
g.write("The estimated cluster number is: " +  str(nClust)+ "\n")
g.close()

# output PCA table with clustering labels
PCA_Data["RealLabel"] = truelabel
PCA_Data["PeakLabel"] = Cluster_Peaks.Label.values
PCA_Data["AllLabel"] = Cluster_All.Label.values
PCA_Data.to_csv(name + "/PCA.csv")

# output normalized gene expression table with clustering labels
S4 = S4.T
S4["RealLabel"] = truelabel
S4["PeakLabel"] = Cluster_Peaks.Label.values
S4["AllLabel"]  = Cluster_All.Label.values
S4.to_csv(name + "/S4.csv")

# output other important results
tSNE.to_csv(name + "/tSNE.csv")

# output clustering results
PH.to_csv(name + "/PersistenceDiagram.csv")
PH_Peaks.to_csv(name + "/PersistenceDiagramPeaks.csv")
Cluster_All.to_csv(name + "/ClusteringLabel.csv")
Cluster_Peaks.to_csv(name + "/ClusteringLabelPeaks.csv")

### Step-11 Further analysis
# intra-cluster stability(pair-wise pearson correlation)
PCA_Data = pd.read_csv(name + "/PCA.csv", header=0, index_col=0, sep=",")
IntraClusterCorrelation(PCA_Data, "RealLabel", name)
IntraClusterCorrelation(PCA_Data, "AllLabel", name)
IntraClusterCorrelation(PCA_Data, "PeakLabel", name)
# downsampling reallabel to the same number of peaklabel
mask = PCA_Data.loc[PCA_Data.PeakLabel == 0,:].shape[0]
sample_ind = list(np.random.choice(list(PCA_Data.index), size=mask, replace=False))
PCA_Data["RealLabel_downsample"] = PCA_Data["RealLabel"]
PCA_Data.loc[sample_ind,"RealLabel_downsample"] = 0
IntraClusterCorrelation(PCA_Data, "RealLabel_downsample", name)

## cluster-specific differentially expressed genes
#S4 = pd.read_csv(name + "/S4.csv", header=0, index_col=0, sep=",")
#tSNE = pd.read_csv(name + "/tSNE.csv", header=0, index_col=0, sep=",")
Markers = ClusterSpecificGeneSet(S4, name, top=10)
