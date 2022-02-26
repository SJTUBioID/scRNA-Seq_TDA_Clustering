# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 15:30:36 2021

@author: cshu
"""

from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import kneighbors_graph
import numpy as np
X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])








starttime = datetime.datetime.now()
KNN = kneighbors_graph(CleanData, 25, mode='connectivity', metric='minkowski', p=2, include_self=False)
KNNmat = KNN.toarray()
endtime = datetime.datetime.now()
print("\nShared Nearest Networt Construction Time:")
print((endtime-starttime))  












starttime = datetime.datetime.now()
Distance = np.round(squareform(pdist(np.array(CleanData), metric='euclidean')),4)
NeighborMatrix = KNNSearch(Distance, k)
endtime = datetime.datetime.now()
print("\nShared Nearest Networt Construction Time:")
print((endtime-starttime))  