#!/usr/bin/env python
# coding: utf-8

# # scRNA-Seq dataset clustering via topologica data analysis (python3)

# ### Step-0 Import packages

# In[1]:



import os
os.chdir("/Users/mac/Desktop/clust/")
from script.functions import *


# ### Step-1 Specify the input dataset
# #### Specify the path of scRNA-Seq dataset. 
# Splatter1_10000 is a synthetic dataset generated from Splatter package, which has 4 equal-sized clusters and 10,000 cells in total.
# 

# In[2]:


file = "Splatter1_10000.txt" 


# #### Specify and create the output path

# In[3]:


name = "Splatter1_10000"


# In[4]:


if not os.path.exists(name):
    os.mkdir(name)


# #### (Optional) Specify the ground-truth cell type label of the dataset. <br/>
# For a synthetic dataset, usually we have a ground-truth cell type label to evaluate the clustering results. <br/>
# However, for a typical true dataset, true cell type label is unknown.

# In[5]:


Label = "Splatter1_10000.label"
Label = pd.read_csv(Label, header=0, index_col=False, sep="\t")
label = Label.values[:,0]
label


# ### Step-2 Data processing
# #### Read dataset into a pandas dataframe
# 

# In[6]:


S1 = pd.read_csv(file, header=0, index_col=0, sep="\s+")
S1.head()


# #### Remove low quality cells and genes
# Remove genes with UMI>0 in less than 3 cells.<br/>
# Remove cells with total UMI<2000 (low-quality).<br/>
# Remove cells with total UMI>10000 (potential doublet).<br/>
# Remove cells with mitochondria UMI contribute >10% of the total transcriptome (low-quality).<br/>

# In[7]:


S2 = Quality_Control(S1)
S2.head()


# #### Normalize Raw UMI count data
# $S3 = np.log(10000*(S2/S2.sum(axis=0))+1) $ <br/>

# In[8]:


S3 = Normalize_Rawdata(S2)
S3.head()


# #### Select genes with highest standard deviation across cells.
# Typically choose 5000~2000 genes

# In[9]:


S4 = Select_HEGenes(S3,ngenes=500)
S4.head()


# #### Scales the expression value of each gene to a N(0,1) distribution.

# In[10]:


S5 = Scaling_Data(S4)
S5.head()


# #### Transpose the scaled matrix to a sample x feature matrix

# In[11]:


S6 = S5.T  
S6.head()


# ### Step-3 Dimension Reduction

# ### Dimension reduction via PCA
# We need to reduce the dimesion S6 matrix by PCA. <br\>
# Choose the first 100 principal components temp temporarily. 

# In[12]:


Dimension = 100
pca = PCA(n_components = Dimension,  random_state = 1) 
PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
PCA_s.head()


# #### Elbowplot displays the variance of each principal component

# In[13]:


elbowplot(S6)


# #### Random permutation helps choose the first n significant principal components
# Permutate each column(gene) of S6 matrix for 100 times, thus we have 100 permuated matricex.
# Perform PCA on each of the 100 permuated matricex and record the variance for PC for each pemutated matrix.
# Then we require that the variance of each PCs of S6 matrix should be higher than the variance of the corresponding PC of permutated matrix.

# In[14]:


pca_permutation(S6, nperm = 100, nPC = 20)


# Based on the elbowplot and pca_permutation test. At least 4 PCs should be kept for this dataset. <br/>
# Users are suggest to try a few times withe number of PCs that is slightly higer than 4.<br/>
# 
# CleanData is the final PCA reduced matrix. <br/>
# Here we choose nDimension = 6 as an example.

# In[15]:


nDimension = 6
CleanData = np.round(PCA_s.iloc[:,0:nDimension],4)
CleanData.head()


# ### Visualization of PCA reduced matrix via tSNE
# PCA reduced matrix may still have tens of dimensions. <br/>
# A further dimesnion reduction for visualization purpose is need to reduce the CleanData matrix into a 2-dimensional tSNE matrix. <br/>
# 

# In[16]:


tsne = TSNE(n_components=2,random_state=1).fit_transform(CleanData)
tsne_embedded = pd.DataFrame(tsne,columns=["tSNE-1","tSNE-2"]) 
tsne_embedded.head()


# ### raw tSNE scatter plot

# In[17]:


simple_tSNE_plot(tsne_embedded)


# ### (optional) tSNE scatter plot with true cell type labels

# In[18]:


tSNE_Plot(tsne_embedded,label)


# ### Step-4 Build neighborhood graph and estimate density for each cell

# ### Labeling each cell with integer ID

# In[19]:


nCell = CleanData.shape[0]
CellID = pd.DataFrame({"Cell":CleanData.index,"ID":np.arange(CleanData.shape[0])})
CellID.head()


# ### Build neighborhood graph and estimate density for each cell

# In[20]:


Edges,Den = SNN_Graph(CleanData,K=30)
Edges.head()
Den.head()


# ### Perform TDA on neighborhood network to find density peaks as representative of each cell type. 

# In[21]:


PH,Cluster,Cluster_U= Persistent_Homology_SNN(CleanData, Edges, Den, name, CellID, tsne_embedded, ind="raw", quiet = True, plot = True)


# In[22]:


PH = PH.sort_values(by="PeakSize" , ascending=False)
print(PH)


# ### (Arbitrary, do not recommend) Remove peaks(clusters) with less than 30 cells 

# In[23]:


if Cluster.shape[1] == 2:
    Cluster[1] = np.ones_like(CellID.ID,dtype="int") * PH.BirthInd[0]
    Cluster["Label"] = Cluster[1]
    PH.PeakSize[0] = CleanData.shape[0]
counts = pd.DataFrame(Cluster.Label.value_counts()) 
cutoff = 30
supress = PH.loc[PH.PeakSize<cutoff,"Barcode"].values
PH_selected = PH.loc[PH.PeakSize>=cutoff,:]
PH_selected.head()


# ### Downsampling to remove small peaks

# In[24]:


# downsampling from 0.9 to 0.1 of original data
for frac in np.linspace(0.9,0.1,num=5):
    print("Downsampling fraction = " + str(frac))
    sample_ind = list(np.random.choice(CellID.ID,int(CellID.shape[0]*frac),replace=False))
    sample_ind.sort()

    subsample  = CleanData.iloc[sample_ind,:]
    subCellID  = CellID.iloc[sample_ind,:]
    subCellID.index = np.arange(int(CellID.shape[0]*frac))  
    sub_tsne_embedded = tsne_embedded.iloc[sample_ind,:]
    sub_tsne_embedded.index = subCellID.index
    subEdges,subDen = SNN_Graph(subsample,K=30)
    subDen.index = subCellID.index
    subPH,subCluster,subCluster_U= Persistent_Homology_SNN(subsample, subEdges, subDen, name, subCellID, sub_tsne_embedded, ind="sampling_"+str(frac), quiet = True, plot = True)
    subPH = subPH.sort_values(by="PeakSize" , ascending=False)
    print(subPH)
    print("***********************")
    print("\n\n\n")


# ### plot TDA clustering predictions after samll peaks removed (high credible points above saddle point). <br/>
# Points below saddle point are unassiged cells. <br/>
# tSNE.Prediction records the clustering results. <br/>
# #### During the downsampling process, we find that there are 4 stable clusters. <br/>
# #### Then the estimated cell type numer K is 4.

# In[25]:


# Since we know there the cluster number is 4
# choose a cutoff to keep the largest 4 peaks in the PH table 
cutoff = 30
supress = PH.loc[PH.PeakSize<cutoff,"Barcode"].values
for c in supress:
    Cluster.loc[Cluster.Label==PH.CurrentCluster[c-1],["Label"]] = -1
Cluster.loc[Cluster.Label ==0,["Label"]] = -1 
  
tsne = pd.DataFrame(tsne_embedded)
tsne["Cell"] = list(CleanData.index)
tsne["Truth"]  =  Label.x
tsne["Prediction"] = Cluster.Label
tSNE_Plot(tsne, tsne.Prediction)  
counts = pd.DataFrame(Cluster.Label.value_counts())  


# ### Replot persistent diagram with cell number for each density peak

# In[26]:


font1 = {'family': 'serif', 'color':  'red', 'weight': 'normal','size': 10, }
idx = PH_selected["PeakSize"].idxmax()
Maximum =  PH_selected.BirthDen.max()
step = Maximum/100
plt.figure(figsize=(8,6),dpi=1200)
plt.plot([0, Maximum], [0, Maximum], ls="--", color="black")
plt.plot(PH_selected.DeathDen,PH_selected.BirthDen, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.2,markersize=6 )
for i in list(PH_selected.index):
    plt.text(PH_selected.DeathDen[i]+step, PH_selected.BirthDen[i]+step, r'${:}$'.format(counts.loc[PH.CurrentCluster[i],"Label"]), fontdict=font1)
plt.title("Persistent Diagram, PeaksNum = " + str(PH_selected .shape[0]))
plt.ylabel("Birth")
plt.xlabel("Death")


# ### Recall rate

# In[27]:


recall  = tsne.loc[tsne.Prediction != -1,:]
recall_rate = np.round(recall.shape[0]/tsne.shape[0],4)
print(recall_rate)


# ### Adjusted Rand Index (optional, when you have true cell type label)

# In[28]:


ARI = np.round(adjusted_rand_score(recall.Truth,recall.Prediction),4)
print(ARI)


# ### Step-5 Reassign the un-assigned cells to each cluster via iterative KNN
# 

# Based on Step-4, a reasonable estimation of Cluster Number K is 4.
# The TDA method explores the topology of a neighborhood graph, assigns points above saddle point to a specific cell type, points below saddle point marked as unassigned. <br/>
# However, when estimating the celluar composition of a given cell population, we need to assign every cell with a cell type labe. This is conducted by an iterative KNN method. <br/> 
# tSNE.Prediction records the clustering results

# In[29]:


tsne = pd.DataFrame(tsne_embedded)
tsne["Cell"] = list(CleanData.index)
tsne["Truth"]  =  Label.x
tsne["Prediction"] = Cluster.Label



Prediction = deepcopy(tsne.Prediction.values)
All_List = list(tsne.index)
UnKnown_List = list(np.where(Prediction<0)[0])
Known_List  = list(np.where(Prediction>0)[0])
    
    
while len(UnKnown_List) >0 :
    UnKnown = CleanData.iloc[UnKnown_List,:]
    Known   = CleanData.iloc[Known_List,:]
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(Known)
    dist, index = nbrs.kneighbors(UnKnown)
    dist  = dist[:,1]
    index = index[:,1]
    known_nearest_index = index[dist.argmin()]
    new_nearest_index = UnKnown_List[dist.argmin()]
    Prediction[new_nearest_index] = Prediction[Known_List[known_nearest_index]] 
    UnKnown_List = list(np.where(Prediction<0)[0])
    Known_List  = list(np.where(Prediction>0)[0])
    #if len(UnKnown_List)%100 == 0:
       # tSNE_Plot(tsne,Prediction) 
tsne["Prediction"] = Prediction
tSNE_Plot(tsne, tsne.Prediction)  
ARI = np.round(adjusted_rand_score(tsne.Truth,tsne.Prediction),4)
print(ARI)


# In[ ]:




