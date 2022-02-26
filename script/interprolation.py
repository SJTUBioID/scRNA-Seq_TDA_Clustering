

nbin = 100
data = tsne
data_max = data.max(axis=0)
data_min = data.min(axis=0)


Grid = []
for i in range(data.shape[1]):
    grid = np.linspace(data_min[i],data_max[i],nbin)
    Grid.append(grid)

GRID = np.array(list(itertools.product(*[Grid[i] for i in range(data.shape[1])])))


from scipy.interpolate import griddata
grid_z0 = griddata(tsne, Den, GRID, method='nearest',fill_value= 0)
grid_z1 = griddata(tsne, Den, GRID, method='linear',fill_value=0)
grid_z2 = griddata(tsne, Den, GRID, method='cubic',fill_value=0)

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8,8)) 
plt.subplot(221)
plt.plot(tsne[:,0], tsne[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.reshape(nbin,nbin).T, origin='lower',cmap="YlGn")
plt.colorbar()
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.reshape(nbin,nbin).T, origin='lower',cmap="YlGn")
plt.colorbar()
plt.title('Linear spline')
plt.subplot(224)
plt.imshow(grid_z2.reshape(nbin,nbin).T, origin='lower',cmap="YlGn")
plt.colorbar()
plt.title('Cubic spline')
plt.gcf().set_size_inches(8, 8)
plt.savefig("tsne_interprolation.png",dpi=600)
plt.close("all")






Den = Den.reshape(10000,1)
Density  = np.hstack((tsne,Den))
scDensity2D1(Density)
scSurface2D1(grid_z0.reshape(nbin,nbin).T)



Density  = np.hstack((GRID,grid_z0.reshape(10000,1)))
Density  = pd.DataFrame(Density,columns=[1,2,"Density"])
PersistentHomologyUnionFindx(Density, name, ind="raw", quiet = True, plot = True)



def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
rng = np.random.default_rng()
points = rng.random((1000, 2))
values = func(points[:,0], points[:,1])


from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')


import matplotlib.pyplot as plt
plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()




nbin = 100
data = points
data_max = data.max(axis=0)
data_min = data.min(axis=0)
Grid = []
for i in range(data.shape[1]):
    grid = np.linspace(data_min[i],data_max[i],nbin)
    Grid.append(grid)

GRID = np.array(list(itertools.product(*[Grid[i] for i in range(data.shape[1])])))


from scipy.interpolate import griddata
grid_z0 = griddata(points, values, GRID, method='nearest')
grid_z1 = griddata(points, values, GRID, method='linear')
grid_z2 = griddata(points, values, GRID, method='cubic')


import matplotlib.pyplot as plt
plt.subplot(221)
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.reshape(nbin,nbin).T, origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z0.reshape(nbin,nbin).T, origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z0.reshape(nbin,nbin).T, origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()




Den = Den.reshape(10000,1)
Density  = np.hstack((tsne,Den))
Density = pd.DataFrame(Density)
scDensity2D1(Density)

values = values.reshape(1000,1)
fuc = np.hstack((points,values))
fuc = pd.DataFrame(fuc)
scDensity2D1(fuc)

