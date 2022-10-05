# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 11:12:07 2020

@author: huchu
"""

import os
os.chdir("C:/Users/huchu/Desktop/DensityCluster/")
from Script.Functions import *

import heapq
import pandas as pd
from scipy import sparse




path = "Jukart_293T/"

class Read_10X_MTX:
    def __init__(self,path):
        if not path.endswith("/"):
            path = path + "/"
        barcode = pd.read_csv(path+"barcodes.tsv", header=None, index_col=None, sep="\t")
        
        if "genes.tsv" in os.listdir(path):
            gene = pd.read_csv(path+"genes.tsv", header=None, index_col=None, sep="\t")
        if "features.tsv" in os.listdir(path):
            gene = pd.read_csv(path+"features.tsv", header=None, index_col=None, sep="\t")
            
        mat = pd.read_csv(path+"matrix.mtx", header=1, index_col=None, sep="\s+",comment='%')
        barcode_names = list(barcode.columns)
        barcode_names[0] = "Barcode"
        barcode.columns = barcode_names

        
        gene_names = list(gene.columns)
        gene_names[0] = "gene_id"
        gene_names[1] = "gene_symbol"
        gene.columns = gene_names

        
        mat.columns = ["gene","barcode","UMI"]
        mat.gene = mat.gene - 1
        mat.barcode = mat.barcode -1
        rows = mat.gene.values
        cols = mat.barcode.values
        v = mat.UMI.values
        sparseM = sparse.csr_matrix((v,(rows,cols)))
    
        self.barcode = barcode
        self.gene = gene
        self.raw = sparseM
        self.cell_id = np.linspace(0,mat.barcode.max(),mat.barcode.max()+1,dtype='int')
        self.gene_id = np.linspace(0,mat.gene.max(),mat.gene.max()+1,dtype='int')


class ReadfromTable:
    def __init__(self,table):
        data = pd.read_csv(table, header=0, index_col=False, sep="\s+")
        
        self.barcode = barcode
        self.gene = gene
        self.raw = mat


class PreProcessing:
    def __init__(self,rawdata):
        self.barcode = rawdata.barcode.copy()
        self.gene = rawdata.gene.copy()
        self.raw = rawdata.raw.copy()
        self.cell_id = rawdata.cell_id.copy()
        self.gene_id = rawdata.gene_id.copy()

        
    def quality_control(self, cell_num = 3, UMI_lowerbound = 2000, UMI_upperbound = 10000, mtPct=0.1):
        raw_count = self.raw.copy()
        gene_symbol = self.gene.gene_symbol.values.copy()
        gene_symbol = gene_symbol[self.gene_id]
        self.gene_symbol = gene_symbol
        
        # Remove cells with high mitochondria expressions
        mt_gene = [True if (s.startswith("MT-") or s.startswith("mt-")) else False  for s in gene_symbol ] 
        self.mtPct = np.array(raw_count[mt_gene,:].sum(axis=0) / raw_count.sum(axis=0))[0,:]
        mt_expr_select = self.mtPct <= mtPct 
        raw_count = raw_count[:,mt_expr_select]
        self.cell_id = self.cell_id[mt_expr_select]
        
        # Select genes
        gene_select = np.array((raw_count>0).sum(axis=1) >= cell_num)[:,0]
        raw_count = raw_count[gene_select,:]
        self.gene_id = self.gene_id[gene_select]
        self.gene_symbol = self.gene_symbol[gene_select]
        # Select cells
        UMI_lower_select = (np.array((raw_count>0).sum(axis = 0) >= UMI_lowerbound)[0,:])
        UMI_upper_select = (np.array((raw_count>0).sum(axis = 0) <= UMI_upperbound)[0,:])
        cell_select =  (UMI_lower_select * UMI_upper_select)
        raw_count = raw_count[:,cell_select]
        self.cell_id = self.cell_id[cell_select]
        
        # To save count matrix after quality control
        self.clean = raw_count        
        
    def normalize(self):
        self.normalize  = sparse.csr_matrix(np.log(10000*(self.clean/self.clean.sum(axis=0))+1))        # Seurat
        #self.normalize = sparse.csr_matrix(np.log2(100000*(self.clean/self.clean.sum(axis=0))+1))      # (SC3 + scTDA) 
 
     
    def select_hvgenes(self,ngene = 500):
        square = self.clean.copy()
        square.data = square.data ** 2
        std = np.array(np.sqrt(square.mean(axis = 1) - np.square(self.clean.mean(axis = 1))))[:,0]
        hvgene_select = heapq.nlargest(ngene,range(std.shape[0]),std.take)
        self.gene_id = self.gene_id[hvgene_select]       
        self.gene_symbol = self.gene_symbol[hvgene_select]
        self.hvmat = pd.DataFrame(self.clean[hvgene_select,:].toarray(), columns = self.cell_id, index = self.gene_id)


name = "Test"
if not os.path.exists(name):
    os.mkdir(name) 
       
raw = Read_10X_MTX(path)
clean = PreProcessing(raw)
clean.quality_control()
clean.normalize()
clean.select_hvgenes(ngene = 500)

S5 = Scaling_Data(clean.hvmat)
S6 = S5.T

### Step-3 Choose the first n Principle Components
Dimension = 100
pca = PCA(n_components=Dimension,  random_state = 1) 
PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
ElbowPlot(S6,name)
PCA_Permutation(S6,name, nperm = 100, nPC = 20)

nDimension = 4
cleanData = np.round(PCA_s.iloc[:,0:nDimension],4)
nCells = cleanData.shape[0]
Density,COOR,bandwidth = sklearnKDE(cleanData)




Runs = 100
starttime = datetime.datetime.now()
test = deepcopy(cleanData)
Result = pd.DataFrame(columns=PersistentDiagramRef.columns)
Result.to_csv(name+"/Permutation_" + str(Runs) + ".txt", sep="\t",header=True,index=False,mode="a")
for i in range(2):
    print(i)
    for j in range(test.shape[1]):
        np.random.shuffle(test.iloc[:,j])   
    density = sklearnKDE_1(test)
    Barcodes, Cluster, Peaks, BarcodeDeath = PersistentHomology(density,name)
    temp = Barcodes.loc[Barcodes.Barcode != 0,:]
    temp.to_csv(name+"/Permutation_" + str(Runs) + ".txt", sep="\t",header=False,index=False,mode="a")

Result =  pd.read_csv(name+"/Permutation_" + str(Runs) + ".txt", header=0, index_col=False, sep="\t")
PeakSizeCutoff = Result.Size.max()
endtime = datetime.datetime.now()
print("\nGaussian Kernel Density Estimation Running Time:")
print((endtime-starttime))  


PersistentDiagramRef["pvalue"] = 1
for i in range(PersistentDiagramRef.shape[0]):
    percentile = stats.percentileofscore(Result.Size, PersistentDiagramRef.Size[i])
    print(100-percentile)
    PersistentDiagramRef.iloc[i,6] = (100-percentile)/100
PersistentDiagramRef.to_csv(name+"/PersistentDiagram_" + str(Runs) + ".txt",sep="\t",header=False,index=False,mode="a")




font1 = {'family': 'serif', 'color':  'skyblue', 'weight': 'normal','size': 10, }
font2 = {'family': 'serif', 'color':  'black', 'weight': 'normal','size': 10, }
idx = Result["Size"].idxmax()
Maximum = Density.Density.max()
Minimum = Density.Density.min()
plt.figure(figsize=(8,6))
plt.plot([0, Maximum], [0, Maximum], ls="--", c=".3", color="black")
plt.plot(Result.Death,Result.Birth, "o", markerfacecolor='red', markeredgecolor='k',markeredgewidth=0.2,markersize=6 )
plt.plot(PersistentDiagramRef.Death,PersistentDiagramRef.Birth, "o", markerfacecolor='skyblue', markeredgecolor='k',markeredgewidth=0.2,markersize=6 )
for i in range(PersistentDiagramRef.shape[0]):
    if PersistentDiagramRef.pvalue[i] <= 0.05:
        plt.text(PersistentDiagramRef.Death[i], PersistentDiagramRef.Birth[i], r'${:.2e}$'.format(PersistentDiagramRef.pvalue[i]), fontdict=font1)
#plt.text(Result.Death[idx],Result.Birth[idx],r'${:.2e}$'.format(Result.Size[idx]), fontdict=font2)
plt.title("Permutated " + str(Runs) + " times")
plt.ylabel("Birth")
plt.xlabel("Death")
plt.savefig(name+"/PersistentDiagramRef_" + str(Runs) + ".txt", dpi=300)
plt.close("all")

