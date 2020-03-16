#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 18:25:22 2019

@author: oem
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import csv
from scipy import signal
from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import h5py

#cd /run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256

#Call and read the file
###############################################################################

#Load the file
filename1='/disk/plasma2/jaa/CB8WAVES/pfd.016000_p000000.h5' 
#Read the file
h5_1st = h5py.File(filename1, 'r')
#Get the fields
###############################################################################
e_field=h5_1st[list(h5_1st.keys())[16]] # Like this the outcome is a group
h_field=h5_1st[list(h5_1st.keys())[17]] # ['hx', 'hy', 'hz']
j_field=h5_1st[list(h5_1st.keys())[18]] # ['jx', 'jy', 'jz']
v_field=h5_1st[list(h5_1st.keys())[24]] # ['vx_e', 'vx_i', 'vy_e', 'vy_i', 'vz_e', 'vz_i']
n_density=h5_1st[list(h5_1st.keys())[22]] #['n_e', 'n_i']
T_tensor=h5_1st[list(h5_1st.keys())[0]] #['Txx_e', 'Txx_i', 'Txy_e', 'Txy_i', 'Txz_e', 'Txz_i', 'Tyy_e', 'Tyy_i', 'Tyz_e', 'Tyz_i', 'Tzz_e', 'Tzz_i']

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


#3D Contour plot
###############################################################################
fig= plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(hx[5,:,:], hy[5,:,:], hz[5,:,:], 50, cmap='viridis')
ax.set_xlabel('hx')
ax.set_ylabel('hy')
ax.set_zlabel('hz')
ax.view_init(60, 35)
fig
#------------------------------------------------------------------------
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(hx[5,:,:], hy[5,:,:], hz[5,:,:], color='black')
ax.set_title('wireframe');
#------------------------------------------------------------------------
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(hx[5,:,:], hy[5,:,:], hz[5,:,:], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('surface');
###############################################################################

fig, ax0 = plt.subplots()
im = plt.pcolor(ex[65,:,:]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('pcolormesh with levels')


#Make operations with the loaded data
###############################################################################

#plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
# ex.shape returns the dimensions of the matrix
mi_me=100;

J_mag = np.sqrt(jx[:,:,:]**2 + jy[:,:,:]**2 + jz[:,:,:]**2)
H_mag = np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2 + hz[:,:,:]**2)
E_mag = np.sqrt(ex[:,:,:]**2 + ey[:,:,:]**2 + ez[:,:,:]**2)
Vi_mag = np.sqrt(vx_i[:,:,:]**2 + vy_i[:,:,:]**2 + vz_i[:,:,:]**2)
Ve_mag = np.sqrt(vx_e[:,:,:]**2 + vy_e[:,:,:]**2 + vz_e[:,:,:]**2)
U_mag =  (n_e[:,:,:]*Ve_mag[:,:,:] + n_i[:,:,:]*Vi_mag[:,:,:]*mi_me)/(n_e[:,:,:] + n_i[:,:,:]*mi_me)


E_par = (ex[:,:,:]*hx[:,:,:] + ey[:,:,:]*hy[:,:,:] + ez[:,:,:]*hz[:,:,:])/H_mag
theta_e = np.arccos((ex[:,:,:]*hx[:,:,:] + ey[:,:,:]*hy[:,:,:] + ez[:,:,:]*hz[:,:,:])/(E_mag*H_mag))
E_perp = E_mag*np.sin(theta_e)

J_par = (jx[:,:,:]*hx[:,:,:] + jy[:,:,:]*hy[:,:,:] + jz[:,:,:]*hz[:,:,:])/H_mag
theta_j = np.arccos((jx[:,:,:]*hx[:,:,:] + jy[:,:,:]*hy[:,:,:] + jz[:,:,:]*hz[:,:,:])/(J_mag*H_mag))
J_perp = J_mag*np.sin(theta_j)

JE_mag=ex[:,:,:]*jx[:,:,:] + ey[:,:,:]*jy[:,:,:] + ez[:,:,:]*jz[:,:,:]



#get the maximum
J_mag_max=np.amax(J_mag[:,:,:])
index_position=np.where(J_mag == J_mag_max) #This returns the index of where that happen

#Get the rms of the current normalise and get the positions
J_mag_rms = np.sqrt(np.mean(J_mag**2))
norm_J_mag=J_mag[:,:,:]/(J_mag_rms)

newArr = J_mag[norm_J_mag > 1]
index_position_2=np.where(norm_J_mag > 1)
plt.pcolor(J_mag[65,:,:])

newArr = J_mag[J_mag > 6*J_mag_rms]
index_position_2=np.where(J_mag > 6*J_mag_rms)



#for i in newArr[:,] 

fig, ax0 = plt.subplots()
im = plt.pcolor(J_mag[:,:,65]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('Current')


#To find the positions where the currents might be

	
# Comparison Operator will be applied to all elements in array
boolArr = J_mag > 6*J_mag_rms
newArr = J_mag[boolArr]

newArr = J_mag[J_mag > 6*J_mag_rms]


#newArr = arr[(arr > 5) & (arr < 20)]


#This is the part made with gorge
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

