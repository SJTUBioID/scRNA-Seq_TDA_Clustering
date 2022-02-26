#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 22:03:35 2021

@author: mac
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:55:42 2021

@author: mac
"""

def PersistentHomologyUnionFindSNN4(CleanData, KNN, Den, name, CellID, umap_embedded, ind="raw", quiet = True, plot = True):
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
    Cluster_U = {}
    
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
        
            
    if plot == True:
        plt.figure(figsize=(8,6))
        plt.plot([0, Den.max()], [0, Den.max()], ls="--", c=".3")
        plt.plot(DeathDensity,BirthDensity, "o", markerfacecolor='coral', markeredgecolor='k',markeredgewidth=0.5,markersize=6 )
        plt.title("Persistent Diagram")
        plt.ylabel("Birth")
        plt.xlabel("Death")
        plt.savefig(name + "/PH_" + str(ind) + ".png",dpi=300)
        plt.close("all")


    ClusterSize[0] = sum(Cluster[1] == Den.argmax())
    PH = pd.DataFrame({"Barcode":Barcode,"BirthInd":BirthIndex,"DeathInd":DeathIndex,"BirthDen":BirthDensity,"DeathDen":DeathDensity,"CurrentCluster":CurrentLevelCluster,"UpperCluster":UpperLeverCluster,"PeakSize":ClusterSize})
    PH["BarcodeLength"] = PH["BirthDen"] - PH["DeathDen"]
    

    U = list(set(PH.UpperCluster))
    U.remove("None")
    U.remove(Den.argmax())
    
    for i in range(len(U)):
        parent_ind = PH.loc[PH.CurrentCluster==U[i],["Barcode"]]
        parent_ind = parent_ind.values.flatten()
        parent_ind = parent_ind[0]
        daughter_ind = PH.loc[PH.UpperCluster==U[i],["Barcode"]]
        daughter_ind = daughter_ind.values.flatten()
        daughter_ind_min = np.array([sum(Cluster_U[i]>0) for i in daughter_ind]).argmin()
        daughter_ind_min = daughter_ind[daughter_ind_min]
        Cluster[parent_ind] = Cluster_U[daughter_ind_min]
    
    
    Cluster = pd.DataFrame(Cluster)
    Cluster["Label"] = Cluster.max(axis=1)   
    endtime = datetime.datetime.now()
    print("\nTopological Data Analysis Running Time:")
    print((endtime-starttime)) 
    return PH,Cluster,Cluster_U