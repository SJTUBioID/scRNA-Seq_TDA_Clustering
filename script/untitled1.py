#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:55:42 2021

@author: mac
"""

def PersistentHomologyUnionFindSNN3(CleanData, KNN, Den, name, CellID, umap_embedded, ind="raw", quiet = True, plot = True):
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
    Edges = deepcopy(KNN)
    
    ID = np.array(CellID.ID.values, dtype="int")
    size = np.ones_like(ID,dtype="int")
    visit = -1*np.ones_like(ID,dtype="int")
    
    BirthIndex = []
    DeathIndex = []
    BirthDensity = []
    DeathDensity = []
    ClusterSize = []
    CurrentLevelCluster = []
    UpperLeverCluster = []
    Barcode=[]
    barcode = 1
    Merge = ["not"]
    Cluster = {1:np.zeros_like(ID,dtype="int")}
   # MergeCluster = {0:visit}
    CurrentLevelCluster.append(Den.argmax())
    UpperLeverCluster.append("None")
    BirthIndex.append(Den.argmax())
    DeathIndex.append(Den.argmin())
    BirthDensity.append(Den.max())
    DeathDensity.append(Den.min())
    ClusterSize.append(0)
    Barcode.append(barcode)
    first = 0

    
    UMAP = pd.concat([umap_embedded,pd.DataFrame(Den)],axis=1)
    
    for i in range(Edges.shape[0]):    
    #for i in range(1400): 
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
                    Label_1 = np.array([root(n) for n in ID ])
                    Label_1[Label_1 != F] = 0
                    Cluster[1] = Label_1
                    print(T,first)
                
                Label = np.array([root(n) for n in ID ])
                Label[visit < 0] = -1    
                scDensity2D(UMAP,name,barcode,Label)
                UMAPPlot(umap_embedded, name, barcode, Label)
                BirthIndex.append(T)
                DeathIndex.append(t)
                BirthDensity.append(Den[T])
                DeathDensity.append(Den[t])
                ClusterSize.append(size[T])
                barcode += 1
                Barcode.append(barcode)
                CleanLabel = deepcopy(Label)
                CleanLabel[CleanLabel != T] = 0
                Cluster[barcode] = CleanLabel
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
                    Label_1 = np.array([root(n) for n in ID ])
                    Label_1[Label_1 != T] = 0
                    Cluster[1] = Label_1
                    print(T,first)               
                
                Label = np.array([root(n) for n in ID ])
                Label[visit < 0] = -1   
                scDensity2D(UMAP,name,barcode,Label)
                UMAPPlot(umap_embedded, name, barcode, Label)
                BirthIndex.append(F)
                DeathIndex.append(t)
                BirthDensity.append(Den[F])
                DeathDensity.append(Den[t])  
                ClusterSize.append(size[F])
                barcode += 1
                Barcode.append(barcode)     
                CleanLabel = deepcopy(Label)
                CleanLabel[CleanLabel != F] = 0
                Cluster[barcode] = CleanLabel
                UpperLeverCluster.append(T)
                CurrentLevelCluster.append(F)
                union(T,F)
      
                
        if Den[F] == Den[T]:
            union(F,T)
        
            
    if plot == True:
        plt.figure(figsize=(8,6))
        plt.plot([0, Den.max()], [0, Den.max()], ls="--", c=".3")
        plt.plot(DeathDensity,BirthDensity, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.5,markersize=6 )
        plt.title("Persistent Diagram")
        plt.ylabel("Birth")
        plt.xlabel("Death")
        plt.savefig(name + "/PH_" + str(ind) + ".png",dpi=300)
        plt.close("all")

    Cluster = pd.DataFrame(Cluster)
    ClusterSize[0] = sum(Cluster[1] == Den.argmax())
    PH = pd.DataFrame({"Barcode":Barcode,"BirthInd":BirthIndex,"DeathInd":DeathIndex,"BirthDen":BirthDensity,"DeathDen":DeathDensity,"CurrentCluster":CurrentLevelCluster,"UpperCluster":UpperLeverCluster,"PeakSize":ClusterSize})
    PH["BarcodeLength"] = PH["BirthDen"] - PH["DeathDen"]
    
    Cluster["Label"] = Cluster.max(axis=1)
    cnts = (Cluster>0).sum(axis=1)
    for i in range(Cluster.shape[0]):
        if cnts[i] >1 :
            temp = Cluster.iloc[i,0:Cluster.shape[1]-1]
            lst = list(temp[temp>0].index)
            current_lst = [PH.CurrentCluster[k-1] for k in lst]
            upper_lst = [PH.UpperCluster[k-1] for k in lst]
            clust_id = [k for k in current_lst if k not in upper_lst]
            Cluster.loc[i,["Label"]] = clust_id[0]
    

    endtime = datetime.datetime.now()
    print("\nTopological Data Analysis Running Time:")
    print((endtime-starttime)) 
    return PH,Cluster