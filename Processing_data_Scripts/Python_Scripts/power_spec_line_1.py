#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 15:21:54 2019

@author: jaa
"""
import numpy as np
import math
from matplotlib import pyplot as plt
import csv
from scipy import signal
from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale

import h5py


#cd /run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256

#Call and read the file
###############################################################################

#Load the file
filename1='/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/1024_stuff/1st/pfd.000600_p000000.h5' 

#Read the file
h5_600_1st = h5py.File(filename1, 'r')
list(h5_600_1st.keys())[16]
"""
['T_1st_single-uid-0x1058210', #0
 'crd[0]-uid-0xd8f670', #1 ['crd[0]']
 'crd[1]-uid-0xd918b0', #2 ['crd[1]']
 'crd[2]-uid-0xd94680', #3 ['crd[2]']
 'crd_nc[0]-uid-0xdaac00', #4 ['crd_nc[0]']
 'crd_nc[1]-uid-0xdaf810',#5 ...
 'crd_nc[2]-uid-0xdb5170', #6 ...
 'crds_gen_x-uid-0xda81f0', #7 []
 'crds_gen_y-uid-0xd95950', #8 []
 'crds_gen_z-uid-0xda9230', #9 []
 'dcrd[0]-uid-0xd97620', #10 []
 'dcrd[1]-uid-0xd9a3f0', #11 []
 'dcrd[2]-uid-0xd9d390', #12 []
 'dcrd_nc[0]-uid-0xda0160', #13 []
 'dcrd_nc[1]-uid-0xda3100', #14 []
 'dcrd_nc[2]-uid-0xda5ed0', #15 []
 'e-uid-0x1a44610', #16
 'h-uid-0x1a3b090', #17
 'j-uid-0x117b840', #18
 'mrc_crds-uid-0xd8d880', #19 []
 'mrc_ddc-uid-0xda71a0', #20 []
 'mrc_domain-uid-0xd6b0f0', #21 []
 'n_1st_single-uid-0xfcfcd0', #22
 'psc-uid-0xf95890', #23 []
 'v_1st_single-uid-0x1a5e6c0'] #24 velocities particles
"""

#Get the fields
###############################################################################
e_field=h5_600_1st[list(h5_600_1st.keys())[16]] # Like this the outcome is a group
h_field=h5_600_1st[list(h5_600_1st.keys())[17]] # ['hx', 'hy', 'hz']
j_field=h5_600_1st[list(h5_600_1st.keys())[18]] # ['jx', 'jy', 'jz']
v_field=h5_600_1st[list(h5_600_1st.keys())[24]] # ['vx_e', 'vx_i', 'vy_e', 'vy_i', 'vz_e', 'vz_i']
n_density=h5_600_1st[list(h5_600_1st.keys())[22]] #['n_e', 'n_i']
T_tensor=h5_600_1st[list(h5_600_1st.keys())[0]] #['Txx_e', 'Txx_i', 'Txy_e', 'Txy_i', 'Txz_e', 'Txz_i', 'Tyy_e', 'Tyy_i', 'Tyz_e', 'Tyz_i', 'Tzz_e', 'Tzz_i']

#coords=h5_600_1st[list(h5_600_1st.keys())[1]]
#list(v_field.keys())


# Electric field components
#-----------------------------------------------------
#Unreference the variables which are in a 3D
ex=e_field['ex']
ex=ex['p0']
ex=ex['3d']

ey=e_field['ey']
ey=ey['p0']
ey=ey['3d']

ez=e_field['ez']
ez=ez['p0']
ez=ez['3d']

# Magnetic field components
#-----------------------------------------------------
hx=h_field['hx']
hx=hx['p0']
hx=hx['3d']

hy=h_field['hy']
hy=hy['p0']
hy=hy['3d']

hz=h_field['hz']
hz=hz['p0']
hz=hz['3d']

# Current field components
#-----------------------------------------------------
jx=j_field['jx']
jx=jx['p0']
jx=jx['3d']

jy=j_field['jy']
jy=jy['p0']
jy=jy['3d']

jz=j_field['jz']
jz=jz['p0']
jz=jz['3d']

# Velocity Ions components
#-----------------------------------------------------
vx_i=v_field['vx_i']
vx_i=vx_i['p0']
vx_i=vx_i['3d']

vy_i=v_field['vy_i']
vy_i=vy_i['p0']
vy_i=vy_i['3d']

vz_i=v_field['vz_i']
vz_i=vz_i['p0']
vz_i=vz_i['3d']

# Velocity electrons components
#-----------------------------------------------------
vx_e=v_field['vx_e']
vx_e=vx_e['p0']
vx_e=vx_e['3d']

vy_e=v_field['vy_e']
vy_e=vy_e['p0']
vy_e=vy_e['3d']

vz_e=v_field['vz_e']
vz_e=vz_e['p0']
vz_e=vz_e['3d']


# Density Ions and electrons
#-----------------------------------------------------
n_i=n_density['n_i']
n_i=n_i['p0']
n_i=n_i['3d']

n_e=n_density['n_e']
n_e=n_e['p0']
n_e=n_e['3d']


# Temperature electron tensor components
#-----------------------------------------------------
Txx_e = T_tensor['Txx_e']
Txx_e=Txx_e['p0']
Txx_e=Txx_e['3d']

Txy_e = T_tensor['Txy_e']
Txy_e=Txy_e['p0']
Txy_e=Txy_e['3d']

Txz_e = T_tensor['Txz_e']
Txz_e=Txz_e['p0']
Txz_e=Txz_e['3d']

Tyy_e = T_tensor['Tyy_e']
Tyy_e=Tyy_e['p0']
Tyy_e=Tyy_e['3d']

Tyz_e = T_tensor['Tyz_e']
Tyz_e=Tyz_e['p0']
Tyz_e=Tyz_e['3d']

Tzz_e = T_tensor['Tzz_e']
Tzz_e=Tzz_e['p0']
Tzz_e=Tzz_e['3d']


# Temperature Ions tensor components
#-----------------------------------------------------
Txx_i = T_tensor['Txx_i']
Txx_i=Txx_i['p0']
Txx_i=Txx_i['3d']

Txy_i = T_tensor['Txy_i']
Txy_i=Txy_i['p0']
Txy_i=Txy_i['3d']

Txz_i = T_tensor['Txz_i']
Txz_i=Txz_i['p0']
Txz_i=Txz_i['3d']

Tyy_i = T_tensor['Tyy_i']
Tyy_i=Tyy_i['p0']
Tyy_i=Tyy_i['3d']

Tyz_i = T_tensor['Tyz_i']
Tyz_i=Tyz_i['p0']
Tyz_i=Tyz_i['3d']

Tzz_i = T_tensor['Tzz_i']
Tzz_i=Tzz_i['p0']
Tzz_i=Tzz_i['3d']
###############################################################################

fig, ax0 = plt.subplots()
im = plt.pcolor(ex[65,:,:]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('pcolormesh with levels')


#plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
Current=np.sqrt(jx[:,:,:]**2 + jy[:,:,:]**2 + jz[:,:,:]**2)
Hfield=np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2 + hz[:,:,:]**2)
Efield=np.sqrt(ex[:,:,:]**2 + ey[:,:,:]**2 + ez[:,:,:]**2)


#get the maximum
Current_max=np.amax(Current[:,:,:])
index_position=np.where(Current == Current_max)

#Get the rms of the current normalise and get the positions
Current_rms = np.sqrt(np.mean(Current**2))
norm_current=Current[:,:,:]/(Current_rms)

newArr = Current[norm_current > 1]
index_position_2=np.where(norm_current > 1)
plt.plot(Current[203,7,:])

newArr = Current[Current > 6*Current_rms]
index_position_2=np.where(Current > 6*Current_rms)




#for i in newArr[:,] 

fig, ax0 = plt.subplots()
im = plt.pcolor(norm_current[500,:,:]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('Current')





#To find the positions where the currents might be

	
# Comparison Operator will be applied to all elements in array
boolArr = Current > 6*Current_rms
newArr = Current[boolArr]

newArr = Current[Current > 6*Current_rms]



#newArr = arr[(arr > 5) & (arr < 20)]

###############################################################################
"""
f5_16 = h5py.File('pfd.00200_p000000.h5', 'r')
list(f5_16.keys())
print("Keys: %s" % f5_16.keys())
f5_2 = h5py.File('pfd.000200_p000000.h5', 'r')
dset_e = f5_2['e-uid-0x1061db0']
ex=dset_e['ex']
p0=ex['p0']
dd_e=p0['3d']
dd_e.shape


plt.pcolor(dd_e[65,:,:]) #plot in the 65 (middle)
plt.pcolor(np.sum(dd_e,axis=0)) #integration over the first dimension

dset_h = f5_16['e-uid-0x1061db0']
ex=dset_e['ex']
p0=ex['p0']
dd_e=p0['3d']
dd_e.shape
plt.pcolor(dd_e[65,:,:]) #plot in the 65 (middle)
plt.pcolor(np.sum(dd_e,axis=0)) #integration over the first dimension
"""
########################################

