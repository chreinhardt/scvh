#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline

# Read data
data = np.loadtxt('scvh_extended_pt_hydrogen_722_helium_278.data', comments='#')

data_logT = data[:,0]
data_logP = data[:,1]
data_logRho = data[:,2]
data_logU = data[:,3]
data_logS = data[:,4]

del data

data_logT_axis = np.sort(np.unique(data_logT))
data_logP_axis = np.sort(np.unique(data_logP))

# Fill into a full grid
logT_array, logP_array = np.meshgrid(data_logT_axis, data_logP_axis, indexing = 'ij')
logRho_array = np.zeros_like(logT_array)
logU_array = np.zeros_like(logT_array)
logS_array = np.zeros_like(logT_array)

for i, logT in enumerate(data_logT_axis):
    for j, logP in enumerate(data_logP_axis):
        found_element_Rho = data_logRho[(data_logT == logT) & (data_logP == logP)]
        logRho_array[i,j] = found_element_Rho[0] if len(found_element_Rho) > 0 else np.nan
        found_element_U = data_logU[(data_logT == logT) & (data_logP == logP)]
        logU_array[i,j] = found_element_U[0] if len(found_element_U) > 0 else np.nan
        found_element_S = data_logS[(data_logT == logT) & (data_logP == logP)]
        logS_array[i,j] = found_element_S[0] if len(found_element_S) > 0 else np.nan

# Extrapolate missing values
for i, logT in enumerate(data_logT_axis):
    filled_indices_Rho = ~np.isnan(logRho_array[i,:])
    interpolatorRho = interp1d(logP_array[i,filled_indices_Rho],logRho_array[i,filled_indices_Rho], kind = 'linear', bounds_error=False, fill_value='extrapolate')
    logRho_array[i,~filled_indices_Rho] = interpolatorRho(logP_array[i,~filled_indices_Rho])
    filled_indices_U = ~np.isnan(logU_array[i,:])
    interpolatorU = interp1d(logP_array[i,filled_indices_U],logU_array[i,filled_indices_U], kind = 'linear', bounds_error=False, fill_value='extrapolate')
    logU_array[i,~filled_indices_U] = interpolatorU(logP_array[i,~filled_indices_U])
    filled_indices_S = ~np.isnan(logS_array[i,:])
    interpolatorS = interp1d(logP_array[i,filled_indices_S],logS_array[i,filled_indices_S], kind = 'linear', bounds_error=False, fill_value='extrapolate')
    logS_array[i,~filled_indices_S] = interpolatorS(logP_array[i,~filled_indices_S])

del interpolatorRho
del interpolatorU
del interpolatorS

interpolatorRho = RectBivariateSpline(data_logT_axis,data_logP_axis,logRho_array,kx = 1, ky = 1, s = 0)
interpolatorU = RectBivariateSpline(data_logT_axis,data_logP_axis,logU_array,kx = 1, ky = 1, s = 0)
interpolatorS = RectBivariateSpline(data_logT_axis,data_logP_axis,logS_array,kx = 1, ky = 1, s = 0)
# interpolated_logT_axis = data_logT_axis
# interpolated_logP_axis = data_logP_axis
interpolated_logT_axis = np.linspace(np.amin(data_logT_axis),np.amax(data_logT_axis),(len(data_logT_axis)-1)*10+1)
interpolated_logP_axis = np.linspace(np.amin(data_logP_axis),np.amax(data_logP_axis),(len(data_logP_axis)-1)*10+1)
interpolated_logT_array, interpolated_logP_array = np.meshgrid(interpolated_logT_axis, interpolated_logP_axis, indexing = 'ij')
interpolated_logRho_array = np.zeros_like(interpolated_logT_array)
interpolated_logU_array = np.zeros_like(interpolated_logT_array)
interpolated_logS_array = np.zeros_like(interpolated_logT_array)

for i, logT in enumerate(interpolated_logT_axis):
    for j, logP in enumerate(interpolated_logP_axis):
        interpolated_logRho_array[i,j] = interpolatorRho(logT,logP)
        interpolated_logU_array[i,j] = interpolatorU(logT,logP)
        interpolated_logS_array[i,j] = interpolatorS(logT,logP)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.plot_surface(interpolated_logT_array,interpolated_logP_array,interpolated_logS_array)
# plt.show()

# Invert into Rho-T tables
min_logRho = np.amin(data_logRho)
max_logRho = np.amax(data_logRho)
inverted_logRho_axis = np.linspace(-21.0,3.0,1001)
inverted_logT_array, inverted_logRho_array = np.meshgrid(interpolated_logT_axis,inverted_logRho_axis, indexing= 'ij')
inverted_logP_array = np.zeros_like(inverted_logT_array)
inverted_logU_array = np.zeros_like(inverted_logT_array)
inverted_logS_array = np.zeros_like(inverted_logT_array)
for i, logT in enumerate(interpolated_logT_axis):
    interpolatorP = interp1d(interpolated_logRho_array[i,:],interpolated_logP_array[i,:], kind = 'linear', bounds_error=False, fill_value='extrapolate')
    inverted_logP_array[i,:] = interpolatorP(inverted_logRho_axis)
    interpolatorU = interp1d(interpolated_logRho_array[i,:],interpolated_logU_array[i,:], kind = 'linear', bounds_error=False, fill_value='extrapolate')
    inverted_logU_array[i,:] = interpolatorU(inverted_logRho_axis)
    interpolatorS = interp1d(interpolated_logRho_array[i,:],interpolated_logS_array[i,:], kind = 'linear', bounds_error=False, fill_value='extrapolate')
    inverted_logS_array[i,:] = interpolatorS(inverted_logRho_axis)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.plot_surface(inverted_logT_array,inverted_logRho_array,inverted_logU_array)
# plt.show()

with open('output.txt','w') as f:
    f.write('# nT = 751 nRho= 1001 (input file: scvh_extended_pt_hydrogen_722_helium_278.data)\n')
    f.write('# logT [K] logRho [g/cc] logP [barye]  logE [erg/g] logS [erg/g/K]\n')
    for i, logT in enumerate(interpolated_logT_axis):
        for j, logRho in enumerate(inverted_logRho_axis):
            f.write('{:.8e} {:.8e} {:.8e} {:.8e} {:.8e}\n'.format(logT,logRho,inverted_logP_array[i,j],inverted_logU_array[i,j],inverted_logS_array[i,j]))

# Debug
print(inverted_logRho_array.shape)

del interpolatorS
interpolatorS = RectBivariateSpline(interpolated_logT_axis,inverted_logRho_axis,inverted_logS_array,kx = 1, ky = 1, s = 0)
fig, ax = plt.subplots(figsize=(16,9))
rho1 = np.log10(202.5599365234375 * 5.6695e-10)
T1 = np.log10(1609.421875)
rho2 = np.log10(3614061.75 * 5.6695e-10)
T2 = np.log10(8450.7529296875)
levels = np.linspace(np.amin(inverted_logS_array),np.amax(inverted_logS_array),101)
levels = np.sort(levels)
ax.contour(inverted_logRho_array,inverted_logT_array,inverted_logS_array,levels)
ax.scatter(rho1,T1)
ax.scatter(rho2,T2)
ax.scatter(data_logRho,data_logT,alpha=0.1)
plt.savefig("S_table.png",dpi=350)

del interpolatorP
interpolatorP = RectBivariateSpline(interpolated_logT_axis,inverted_logRho_axis,inverted_logP_array,kx = 1, ky = 1, s = 0)
fig, ax = plt.subplots(figsize=(16,9))
rho1 = np.log10(202.5599365234375 * 5.6695e-10)
T1 = np.log10(1609.421875)
rho2 = np.log10(3614061.75 * 5.6695e-10)
T2 = np.log10(8450.7529296875)
levels = np.linspace(np.amin(inverted_logP_array),np.amax(inverted_logP_array),101)
levels = np.sort(levels)
ax.contour(inverted_logRho_array,inverted_logT_array,inverted_logP_array,levels)
ax.scatter(rho1,T1)
ax.scatter(rho2,T2)
ax.scatter(data_logRho,data_logT,alpha=0.1)
plt.savefig("P_table.png",dpi=350)


del interpolatorU
interpolatorU = RectBivariateSpline(interpolated_logT_axis,inverted_logRho_axis,inverted_logU_array,kx = 1, ky = 1, s = 0)
fig, ax = plt.subplots(figsize=(16,9))
rho1 = np.log10(202.5599365234375 * 5.6695e-10)
T1 = np.log10(1609.421875)
rho2 = np.log10(3614061.75 * 5.6695e-10)
T2 = np.log10(8450.7529296875)
levels = np.linspace(np.amin(inverted_logU_array),np.amax(inverted_logU_array),101)
levels = np.sort(levels)
ax.contour(inverted_logRho_array,inverted_logT_array,inverted_logU_array,levels)
ax.scatter(rho1,T1)
ax.scatter(rho2,T2)
ax.scatter(data_logRho,data_logT,alpha=0.1)
plt.savefig("U_table.png",dpi=350)
plt.show()


