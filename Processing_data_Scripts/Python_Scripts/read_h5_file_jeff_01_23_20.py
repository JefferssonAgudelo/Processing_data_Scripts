#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 11:05:58 2020

This is a script that can be used to open a single .h5 file. 

@author: jaa
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import glob, os
import gc

gc.enable()


###############################################################################
# Load the file

path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1'  # Put your path here

file=sorted(glob.glob(os.path.join(path, '*.h5')))
a=1


filename1= file[a]
print(filename1)
h5_1st = h5py.File(filename1, 'r')

e_field=h5_1st[list(h5_1st.keys())[16]] # Like this the outcome is a group
h_field=h5_1st[list(h5_1st.keys())[17]] # ['hx', 'hy', 'hz']
j_field=h5_1st[list(h5_1st.keys())[18]] # ['jx', 'jy', 'jz']
v_field=h5_1st[list(h5_1st.keys())[24]] # ['vx_e', 'vx_i', 'vy_e', 'vy_i', 'vz_e', 'vz_i']
n_density=h5_1st[list(h5_1st.keys())[22]] #['n_e', 'n_i']
T_tensor=h5_1st[list(h5_1st.keys())[0]] #['Txx_e', 'Txx_i', 'Txy_e', 'Txy_i', 'Txz_e', 'Txz_i', 'Tyy_e', 'Tyy_i', 'Tyz_e', 'Tyz_i', 'Tzz_e', 'Tzz_i']

ex=e_field['ex'] ; ex=ex['p0'] ; ex=ex['3d'] ; #ex=ex[:,:,:] ; ex=np.transpose(ex, (1, 2, 0))
ey=e_field['ey'] ; ey=ey['p0'] ; ey=ey['3d'] ; #ey=ey[:,:,:] ; ey=np.transpose(ey, (1, 2, 0))
ez=e_field['ez'] ; ez=ez['p0'] ; ez=ez['3d'] ; #ez=ez[:,:,:] ; ez=np.transpose(ez, (1, 2, 0))

# If the file is too heavy this will take some time. It is advisable to load just some of the variables, for example just the magnetic field. 
# In that case, comment the other variables

# Magnetic field components
#-----------------------------------------------------
hx=h_field['hx'] ; hx=hx['p0'] ; hx=hx['3d'] ; #hx=hx[:,:,:] ; hx=np.transpose(hx, (1, 2, 0))
hy=h_field['hy'] ; hy=hy['p0'] ; hy=hy['3d'] ; #hy=hy[:,:,:] ; hy=np.transpose(hy, (1, 2, 0))
hz=h_field['hz'] ; hz=hz['p0'] ; hz=hz['3d'] ; #hz=hz[:,:,:] ; hz=np.transpose(hz, (1, 2, 0))

# Current field components
#-----------------------------------------------------
jx=j_field['jx'] ; jx=jx['p0'] ; jx=jx['3d'] ; #jx=jx[:,:,:] ; jx=np.transpose(jx, (1, 2, 0))
jy=j_field['jy'] ; jy=jy['p0'] ; jy=jy['3d'] ; #jy=jy[:,:,:] ; jy=np.transpose(jy, (1, 2, 0))
jz=j_field['jz'] ; jz=jz['p0'] ; jz=jz['3d'] ; #jz=jz[:,:,:] ; jz=np.transpose(jz, (1, 2, 0))

# Velocity Ions components
#-----------------------------------------------------
vx_i=v_field['vx_i'] ; vx_i=vx_i['p0'] ; vx_i=vx_i['3d'] #; vx_i=vx_i[:,:,:] ; vx_i=np.transpose(vx_i, (1, 2, 0))
vy_i=v_field['vy_i'] ; vy_i=vy_i['p0'] ; vy_i=vy_i['3d'] #; vy_i=vy_i[:,:,:] ; vy_i=np.transpose(vy_i, (1, 2, 0))
vz_i=v_field['vz_i'] ; vz_i=vz_i['p0'] ; vz_i=vz_i['3d'] #; vz_i=vz_i[:,:,:] ; vz_i=np.transpose(vz_i, (1, 2, 0))

# Velocity electrons components
#-----------------------------------------------------
vx_e=v_field['vx_e'] ; vx_e=vx_e['p0'] ; vx_e=vx_e['3d'] #; vx_e=vx_e[:,:,:] ; vx_e=np.transpose(vx_e, (1, 2, 0))
vy_e=v_field['vy_e'] ; vy_e=vy_e['p0'] ; vy_e=vy_e['3d'] #; vy_e=vy_e[:,:,:] ; vy_e=np.transpose(vy_e, (1, 2, 0))
vz_e=v_field['vz_e'] ; vz_e=vz_e['p0'] ; vz_e=vz_e['3d'] #; vz_e=vz_e[:,:,:] ; vz_e=np.transpose(vz_e, (1, 2, 0))

# Density Ions and electrons
#-----------------------------------------------------
n_i=n_density['n_i'] ; n_i=n_i['p0'] ; n_i=n_i['3d'] ; #n_i=n_i[:,:,:] ; n_i=np.transpose(n_i, (1, 2, 0))
n_e=n_density['n_e'] ; n_e=n_e['p0'] ; n_e=n_e['3d'] ; #n_e=n_e[:,:,:] ; n_e=np.transpose(n_e, (1, 2, 0))

# Temperature electron tensor components T_1st_single
#-----------------------------------------------------
Txx_e = T_tensor['Txx_e'] ; Txx_e=Txx_e['p0'] ; Txx_e=Txx_e['3d'] ; #Txx_e=np.transpose(Txx_e, (1, 2, 0))
Txy_e = T_tensor['Txy_e'] ; Txy_e=Txy_e['p0'] ; Txy_e=Txy_e['3d'] ; #Txy_e=np.transpose(Txy_e, (1, 2, 0))
Txz_e = T_tensor['Txz_e'] ; Txz_e=Txz_e['p0'] ; Txz_e=Txz_e['3d'] ; #Txz_e=np.transpose(Txz_e, (1, 2, 0))
Tyy_e = T_tensor['Tyy_e'] ; Tyy_e=Tyy_e['p0'] ; Tyy_e=Tyy_e['3d'] ; #Tyy_e=np.transpose(Tyy_e, (1, 2, 0))
Tyz_e = T_tensor['Tyz_e'] ; Tyz_e=Tyz_e['p0'] ; Tyz_e=Tyz_e['3d'] ; #Tyz_e=np.transpose(Tyz_e, (1, 2, 0))
Tzz_e = T_tensor['Tzz_e'] ; Tzz_e=Tzz_e['p0'] ; Tzz_e=Tzz_e['3d'] ; #Tzz_e=np.transpose(Tzz_e, (1, 2, 0))
# Temperature Ions tensor components
#-----------------------------------------------------
Txx_i = T_tensor['Txx_i'] ; Txx_i=Txx_i['p0'] ; Txx_i=Txx_i['3d'] ; #Txx_i=np.transpose(Txx_i, (1, 2, 0))
Txy_i = T_tensor['Txy_i'] ; Txy_i=Txy_i['p0'] ; Txy_i=Txy_i['3d'] ; #Txy_i=np.transpose(Txy_i, (1, 2, 0))
Txz_i = T_tensor['Txz_i'] ; Txz_i=Txz_i['p0'] ; Txz_i=Txz_i['3d'] ; #Txz_i=np.transpose(Txz_i, (1, 2, 0))
Tyy_i = T_tensor['Tyy_i'] ; Tyy_i=Tyy_i['p0'] ; Tyy_i=Tyy_i['3d'] ; #Tyy_i=np.transpose(Tyy_i, (1, 2, 0))
Tyz_i = T_tensor['Tyz_i'] ; Tyz_i=Tyz_i['p0'] ; Tyz_i=Tyz_i['3d'] ; #Tyz_i=np.transpose(Tyz_i, (1, 2, 0))
Tzz_i = T_tensor['Tzz_i'] ; Tzz_i=Tzz_i['p0'] ; Tzz_i=Tzz_i['3d'] ; #Tzz_i=np.transpose(Tzz_i, (1, 2, 0))
###############################################################################


############################################################################################################
mi_me=100;
H_mag = np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2 + hz[:,:,:]**2)
J_mag = np.sqrt(jx[:,:,:]**2 + jy[:,:,:]**2 + jz[:,:,:]**2)
E_mag = np.sqrt(ex[:,:,:]**2 + ey[:,:,:]**2 + ez[:,:,:]**2)
Vi_mag = np.sqrt(vx_i[:,:,:]**2 + vy_i[:,:,:]**2 + vz_i[:,:,:]**2)
Ve_mag = np.sqrt(vx_e[:,:,:]**2 + vy_e[:,:,:]**2 + vz_e[:,:,:]**2)
U_mag =  (n_e[:,:,:]*Ve_mag + n_i[:,:,:]*Vi_mag*mi_me)/(n_e[:,:,:] + n_i[:,:,:]*mi_me)


# Temperatures must be independet of the mass and they are
#------------------------------------------------------------------------------
Ti = (Txx_i[:,:,:] + Tyy_i[:,:,:] + Tzz_i[:,:,:])/3
Te = (Txx_e[:,:,:] + Tyy_e[:,:,:] + Tzz_e[:,:,:])/3
Vith = np.sqrt(2)*np.sqrt(np.abs(Ti))
Veth = np.sqrt(2*mi_me)*np.sqrt(np.abs(Te))

Ti_R = ((Txx_i[:,:,:] - vx_i[:,:,:]**2) + (Tyy_i[:,:,:] - vy_i[:,:,:]**2) + (Tzz_i[:,:,:]- vz_i[:,:,:]**2))/3
Te_R = ((Txx_e[:,:,:] - vx_e[:,:,:]**2) + (Tyy_e[:,:,:] - vy_e[:,:,:]**2) + (Tzz_e[:,:,:]- vz_e[:,:,:]**2))/3
Vith_R = np.sqrt(2)*np.sqrt(np.abs(Ti_R))
Veth_R = np.sqrt(2*mi_me)*np.sqrt(np.abs(Te_R))

#------------------------------------------------------------------------------


####################################################################################################################
gc.collect()



# Make some plots
################################################################################################
fig, ax0 = plt.subplots()
im = plt.pcolor(H_mag[:,:,0]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('pcolormesh with levels')        

fig, ax0 = plt.subplots()
im = plt.pcolor(J_mag[:,:,0]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('pcolormesh with levels')        

fig, ax0 = plt.subplots()
im = plt.pcolor(E_mag[:,:,0]) #plot in the 65 (middle)
fig.colorbar(im, ax=ax0)
ax0.set_title('E_mag')    

#################################################################################################################





















