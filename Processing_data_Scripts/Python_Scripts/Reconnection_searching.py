#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:56:25 2019

@author: jaa
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import glob, os
import gc
import pylab

gc.enable()

import pprint


import resource

from mpl_toolkits import mplot3d
from scipy import stats
#from scipy.stats import norm
#from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
#import csv
#import math

#from __future__ import print_function
resource.getrusage(resource.RUSAGE_SELF).ru_maxrss



###############################################################################
# Load the file
#path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'
path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1'

path = '/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one'

file=sorted(glob.glob(os.path.join(path, '*.h5')))
a=5 #choose the file
a=1
#Utotal_av=[]; U_elec_av=[]; U_mag_av =[]; Ki_av =[]; Ke_av =[]; Ki_th_av =[]; Ke_th_av=[]
#J_rms_T=[]; B_rms_T=[]; Vi_rms_T=[]

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

# Magnetic field components
#-----------------------------------------------------
hx=h_field['hx'] ; hx=hx['p0'] ; hx=hx['3d'] ; hx=hx[:,:,:] ; hx=np.transpose(hx, (1, 2, 0))
hy=h_field['hy'] ; hy=hy['p0'] ; hy=hy['3d'] ; hy=hy[:,:,:] ; hy=np.transpose(hy, (1, 2, 0))
hz=h_field['hz'] ; hz=hz['p0'] ; hz=hz['3d'] ; hz=hz[:,:,:] ; hz=np.transpose(hz, (1, 2, 0))

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


#####################################################################
# This is how to pass from np.array to pd.dataframe
data = np.array([['','Col1','Col2'],['Row1',1,2],['Row2',3,4]])
pd.DataFrame(data=data[1:,1:],    # values
index=data[1:,0],    # 1st column as index
columns=data[0,1:])  # 1st row as the column names


############################################################################################################
mi_me=100;
H_mag = np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2 + hz[:,:,:]**2)
J_mag = np.sqrt(jx[:,:,:]**2 + jy[:,:,:]**2 + jz[:,:,:]**2)
E_mag = np.sqrt(ex[:,:,:]**2 + ey[:,:,:]**2 + ez[:,:,:]**2)
Vi_mag = np.sqrt(vx_i[:,:,:]**2 + vy_i[:,:,:]**2 + vz_i[:,:,:]**2)
Ve_mag = np.sqrt(vx_e[:,:,:]**2 + vy_e[:,:,:]**2 + vz_e[:,:,:]**2)
U_mag =  (n_e[:,:,:]*Ve_mag + n_i[:,:,:]*Vi_mag*mi_me)/(n_e[:,:,:] + n_i[:,:,:]*mi_me)

# Temperature must be independet of the mass and they are
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

E_par = (ex[:,:,:]*hx[:,:,:] + ey[:,:,:]*hy[:,:,:] + ez[:,:,:]*hz[:,:,:])/H_mag
theta_e = np.arccos((ex[:,:,:]*hx[:,:,:] + ey[:,:,:]*hy[:,:,:] + ez[:,:,:]*hz[:,:,:])/(E_mag*H_mag))
E_perp = E_mag*np.sin(theta_e)
J_par = (jx[:,:,:]*hx[:,:,:] + jy[:,:,:]*hy[:,:,:] + jz[:,:,:]*hz[:,:,:])/H_mag
theta_j = np.arccos((jx[:,:,:]*hx[:,:,:] + jy[:,:,:]*hy[:,:,:] + jz[:,:,:]*hz[:,:,:])/(J_mag*H_mag))
J_perp = J_mag*np.sin(theta_j)
JE_mag=ex[:,:,:]*jx[:,:,:] + ey[:,:,:]*jy[:,:,:] + ez[:,:,:]*jz[:,:,:]


Vi_m_ExB2_x = (vx_i[:,:,:]-(ey[:,:,:]*hz[:,:,:]-ez[:,:,:]*hy[:,:,:])/H_mag**2)
Vi_m_ExB2_y = (vy_i[:,:,:]-(ez[:,:,:]*hx[:,:,:]-ex[:,:,:]*hz[:,:,:])/H_mag**2)
Vi_m_ExB2_z = (vz_i[:,:,:]-(ex[:,:,:]*hy[:,:,:]-ey[:,:,:]*hx[:,:,:])/H_mag**2)
Vi_m_ExB2 = np.sqrt(Vi_m_ExB2_x*Vi_m_ExB2_x + Vi_m_ExB2_y*Vi_m_ExB2_y + Vi_m_ExB2_z*Vi_m_ExB2_z)
del(Vi_m_ExB2_x) ; del(Vi_m_ExB2_y) ; del(Vi_m_ExB2_z)  


Ve_m_ExB2_x = (vx_e[:,:,:]-(ey[:,:,:]*hz[:,:,:]-ez[:,:,:]*hy[:,:,:])/H_mag**2) 
Ve_m_ExB2_y = (vy_e[:,:,:]-(ez[:,:,:]*hx[:,:,:]-ex[:,:,:]*hz[:,:,:])/H_mag**2) 
Ve_m_ExB2_z = (vz_e[:,:,:]-(ex[:,:,:]*hy[:,:,:]-ey[:,:,:]*hx[:,:,:])/H_mag**2)
Ve_m_ExB2 = np.sqrt(Ve_m_ExB2_x*Ve_m_ExB2_x + Ve_m_ExB2_y*Ve_m_ExB2_y + Ve_m_ExB2_z*Ve_m_ExB2_z)
del(Ve_m_ExB2_x) ; del(Ve_m_ExB2_y) ; del(Ve_m_ExB2_z)

U_m_ExB_B2_x = ((mi_me*n_i[:,:,:]*vx_i[:,:,:]+n_e[:,:,:]*vx_e[:,:,:])/(n_e[:,:,:] + mi_me*n_i[:,:,:]) - (ey[:,:,:]*hz[:,:,:]-ez[:,:,:]*hy[:,:,:])/H_mag**2) 
U_m_ExB_B2_y = ((mi_me*n_i[:,:,:]*vy_i[:,:,:]+n_e[:,:,:]*vy_e[:,:,:])/(n_e[:,:,:] + mi_me*n_i[:,:,:]) - (ez[:,:,:]*hx[:,:,:]-ex[:,:,:]*hz[:,:,:])/H_mag**2)
U_m_ExB_B2_z = ((mi_me*n_i[:,:,:]*vz_i[:,:,:]+n_e[:,:,:]*vz_e[:,:,:])/(n_e[:,:,:] + mi_me*n_i[:,:,:]) - (ex[:,:,:]*hy[:,:,:]-ey[:,:,:]*hx[:,:,:])/H_mag**2)
U_m_ExB_B2 = np.sqrt(U_m_ExB_B2_x*U_m_ExB_B2_x + U_m_ExB_B2_y*U_m_ExB_B2_y + U_m_ExB_B2_z*U_m_ExB_B2_z)
del(U_m_ExB_B2_x) ; del(U_m_ExB_B2_y) ; del(U_m_ExB_B2_z)



####################################################################################################################
gc.collect()




nx=len(hx[:,0,0]) # This takes a lot of time better just put by hand
ny=len(hx[0,:,0])
nz=len(hx[0,0,:])

#Divergence
###############################################################################
divB=np.float32(np.zeros((nx-1,ny-1,nz-1)))
divBx=np.float32(np.zeros((nx-1,ny-1,nz-1)))
divBy=np.float32(np.zeros((nx-1,ny-1,nz-1)))
divBz=np.float32(np.zeros((nx-1,ny-1,nz-1)))

divB2=np.float32(np.zeros((nx-1,ny-1,nz-1)))
divBx2=np.float32(np.zeros((nx-1,ny-1,nz-1)))
divBy2=np.float32(np.zeros((nx-1,ny-1,nz-1)))
divBz2=np.float32(np.zeros((nx-1,ny-1,nz-1)))

# Divergence in the usal way.-------------------------------------------------
for k in range(nz-1):
    for j in range(ny-1):
        for i in range(nx-1):
            divBx[i,j,k]=(hx[i+1,j,k]-hx[i,j,k])
            divBy[i,j,k]=(hy[i,j+1,k]-hy[i,j,k])
            divBz[i,j,k]=(hz[i,j,k+1]-hz[i,j,k])
divB = divBx + divBy + divBz


#Divergence in the stagered yee grid------------------------------------------
for k in range(nz-2):
    for j in range(ny-2):
        for i in range(nx-2):
            divBx2[i+1,j+1,k+1]=(hx[i+2,j+1,k+1]-hx[i+1,j+1,k+1])
            divBy2[i+1,j+1,k+1]=(hy[i+1,j+2,k+1]-hy[i+1,j+1,k+1])
            divBz2[i+1,j+1,k+1]=(hz[i+1,j+1,k+2]-hz[i+1,j+1,k+1])
divB2 = divBx2 + divBy2 + divBz2

    
hx[i+1,j,k]-hx[i-1,j,k]
hx[1,0,0]-hx[-1,0,0]



###############################################################################

curlBx=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlBy=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlBz=np.float32(np.zeros((nx-1,ny-1,nz-1)))
magCurlB_2=np.float32(np.zeros((nx-1,ny-1,nz-1)))

curlVix=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlViy=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlViz=np.float32(np.zeros((nx-1,ny-1,nz-1)))
magCurlVi_2=np.float32(np.zeros((nx-1,ny-1,nz-1)))

curlVex=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlVey=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlVez=np.float32(np.zeros((nx-1,ny-1,nz-1)))
magCurlVe_2=np.float32(np.zeros((nx-1,ny-1,nz-1)))



mi_me=100; #Should it be with N=1?
U_x =  (vx_e + mi_me*vx_i)/(1+mi_me)
U_y =  (vy_e + mi_me*vy_i)/(1+mi_me)
U_z =  (vz_e + mi_me*vz_i)/(1+mi_me)  

curlUx=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlUy=np.float32(np.zeros((nx-1,ny-1,nz-1)))
curlUz=np.float32(np.zeros((nx-1,ny-1,nz-1)))
magCurlU_2=np.float32(np.zeros((nx-1,ny-1,nz-1)))


#distance = [[[0 for k in range(nz)] for j in range(ny)] for i in range(nx)]
#pprint.pprint(distance)


# Calculate the curl of the magnetic field
for k in range(nz-1):
    for j in range(ny-1):
        for i in range(nx-1):
            #curlBx[i,j,k]=(hz[i,j+1,k]-hz[i,j,k]) - (hy[i,j,k+1]-hy[i,j,k])
            #curlBy[i,j,k]=(hx[i,j,k+1]-hx[i,j,k]) - (hz[i+1,j,k]-hz[i,j,k])
            curlBz[i,j,k]=(hy[i+1,j,k]-hy[i,j,k]) - (hx[i,j+1,k]-hx[i,j,k])
#magCurlB_2 = curlBx**2 + curlBy**2 + curlBz**2 # Joule heating

for k in range(nz-1):
    for j in range(ny-1):
        for i in range(nx-1):
            #curlVix[i,j,k]=(vz_i[i,j+1,k]-vz_i[i,j,k]) - (vy_i[i,j,k+1]-vy_i[i,j,k])
            #curlViy[i,j,k]=(vx_i[i,j,k+1]-vx_i[i,j,k]) - (vz_i[i+1,j,k]-vz_i[i,j,k])
            curlViz[i,j,k]=(vy_i[i+1,j,k]-vy_i[i,j,k]) - (vx_i[i,j+1,k]-vx_i[i,j,k])
#magCurlvi_2 = #curlVix**2 + curlViy**2 + curlViz**2

for k in range(nz-1):
    for j in range(ny-1):
        for i in range(nx-1):
            #curlVex[i,j,k]=(vz_e[i,j+1,k]-vz_e[i,j,k]) - (vy_e[i,j,k+1]-vy_e[i,j,k])
            #curlVey[i,j,k]=(vx_e[i,j,k+1]-vx_e[i,j,k]) - (vz_e[i+1,j,k]-vz_e[i,j,k])
            curlVez[i,j,k]=(vy_e[i+1,j,k]-vy_e[i,j,k]) - (vx_e[i,j+1,k]-vx_e[i,j,k])
#magCurlve_2 = #curlVex**2 + curlVey**2 + curlVez**2


for k in range(nz-1):
    for j in range(ny-1):
        for i in range(nx-1):
            curlUx[i,j,k]=(U_z[i,j+1,k]-U_z[i,j,k]) - (U_y[i,j,k+1]-U_y[i,j,k])
            #curlUy[i,j,k]=(U_x[i,j,k+1]-U_x[i,j,k]) - (U_z[i+1,j,k]-U_z[i,j,k])
            #curlUz[i,j,k]=(U_y[i+1,j,k]-U_y[i,j,k]) - (U_x[i,j+1,k]-U_x[i,j,k])
#magCurlU_2 = #curlUx**2 + curlUy**2 + curlUz**2


#To check where where the vorticity is opposite 
            
curlBzcurlViz = curlBz*curlViz             
signBVi=np.sign(curlBzcurlViz)
curlBzcurlVez = curlBz*curlVez   
signBVe=np.sign(curlBzcurlVez)
################################################################################################

# Get the directory
os.getcwd()            
os.chdir('/disk/plasma2/jaa/CB8WAVES/CB8waves_04/Current_along_box')

#Plots of the curls
################################################################################################
for i in range(10):
    nn= int(i)
    fig, ax0 = plt.subplots()
    im = plt.pcolor(J_sheet_2_Uri[:,:,nn]) #plot in the 65 (middle)
    fig.colorbar(im, ax=ax0)
    ax0.set_title('Current_sheet')
    fig.savefig("current"+ str(nn)+".png", bbox_inches='tight')


###############################################################################
pylab.pcolor(divB[:,:,15], cmap='bwr')
pylab.colorbar()
pylab.title('Div B')
#pylab.savefig('120_curlBzVi.png')
pylab.show() 

pylab.pcolor(divB2[:,:,15], cmap='bwr')
pylab.colorbar()
pylab.title('Div B')
#pylab.savefig('120_curlBzVi.png')
pylab.show() 

pylab.pcolor(divBx[:,:,100], cmap='bwr')
pylab.colorbar()
pylab.title('Div B')
#pylab.savefig('120_curlBzVi.png')
pylab.show() 
###############################################################################


pylab.pcolor(signBVi[:,:,0], cmap='bwr')
pylab.colorbar()
pylab.title('Curl Vi')
#pylab.savefig('120_curlBzVi.png')
pylab.show() 

f2=pylab.pcolor(signBVe[:,:,0], cmap='bwr')
pylab.colorbar()
pylab.title('Curl Ve')
#pylab.savefig('120_curlBzVe.png')
pylab.show() 

f3=pylab.pcolor(curlVez[:,:,0], cmap='bwr')
pylab.colorbar()
pylab.title('Curl Vez')
pylab.savefig('120_curlVez.png')
pylab.show() 



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

# How to save and load the data save as npy
#################################################################################################################

np.savez('curls120.npz', curlBz=curlBz, curlViz=curlViz, curlVez=curlVez, signBVi=signBVi, signBVe=signBVe, llow_pickle=False)
data=np.load('curls120.npz')

curlBz1 = data['curlBz']
curlViz1 = data['curlViz']
curlVez1 = data['curlVez']
signBVi1 = data['signBVi']
signBVe1 = data['signBVe']
#################################################################################################################


#####################################################################################################
# The current given by the simulation and the one calculate with the curl of B dony exactly
# match. 
# To save the outputs and open them in paraview I need to work with the coordinates in the origina file
# otherwise there are no volume values in paraview 
#######################################################################################################
from Xdmf import *
# Create a Domain
domain = XdmfDomain.New()
hf = h5py.File('curl_B_U.h5', 'w')
hf.create_dataset('curlBx', data=curlBx)
hf.create_dataset('curlBy', data=curlBy)
hf.create_dataset('curlBz', data=curlBz)
hf.create_dataset('curlUx', data=curlUx)
hf.create_dataset('curlUy', data=curlUy)
hf.create_dataset('curlUz', data=curlUz)
###########################################################################################################



#Calculating statistics
###############################################################################################################

H_mag_mean = np.mean(H_mag); H_mag_std = np.std(H_mag);
E_mag_mean = np.mean(E_mag); E_mag_std = np.std(E_mag);
J_mag_mean = np.mean(J_mag); J_mag_std = np.std(J_mag);
Vi_mag_mean = np.mean(Vi_mag); Vi_mag_std = np.std(Vi_mag);
Ve_mag_mean = np.mean(Ve_mag); Ve_mag_std = np.std(Ve_mag);
U_mag_mean = np.mean(U_mag); U_mag_std = np.std(U_mag);

Ti_mean = np.mean(Ti); Ti_std = np.std(Ti);
Te_mean = np.mean(Te); Te_std = np.std(Te);
Vith_mean = np.mean(Vith); Vith_std = np.std(Vith);
Veth_mean = np.mean(Veth); Veth_std = np.std(Veth);

Ti_mean_R = np.mean(Ti_R); Ti_std_R = np.std(Ti_R);
Te_mean_R = np.mean(Te_R); Te_std_R = np.std(Te_R);
Vith_mean_R = np.mean(Vith_R); Vith_std_R = np.std(Vith_R);
Veth_mean_R = np.mean(Veth_R); Veth_std_R = np.std(Veth_R);

E_par_mean = np.mean(E_par); E_par_std = np.std(E_par);
E_perp_mean = np.mean(E_perp); E_perp_std = np.std(E_perp);
J_par_mean = np.mean(J_par); J_par_std = np.std(J_par);
J_perp_mean = np.mean(J_perp); J_perp_std = np.std(J_perp);


# Calculating the thrseholds for the isosurfaces an representatives values
Hmag_threshold = H_mag_mean + 4*H_mag_std
Emag_threshold = E_mag_mean + 4*E_mag_std
Jmag_threshold = J_mag_mean + 4*J_mag_std
Vi_threshold = Vi_mag_mean + 4*Vi_mag_std
Ve_threshold = Ve_mag_mean + 4*Ve_mag_std
U_threshold = U_mag_mean + 4*U_mag_std

Ti_threshold = Ti_mean + 4*Ti_std
Te_threshold = Te_mean + 4*Te_std
Vith_threshold = Vith_mean + 4*Vith_std
Veth_threshold = Veth_mean + 4*Veth_std
Ti_threshold_R = Ti_mean_R + 4*Ti_std_R
Te_threshold_R = Te_mean_R + 4*Te_std_R
Vith_threshold_R = Vith_mean_R + 4*Vith_std_R
Veth_threshold_R = Veth_mean_R + 4*Veth_std_R

Epar_threshold = E_par_mean + 4*E_par_std
Eperp_threshold = E_perp_mean + 4*E_perp_std
Jpar_threshold = J_par_mean + 4*J_par_std
Jperp_threshold = J_perp_mean + 4*J_perp_std

         
os.chdir('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/')
#os.chdir('/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one')
# Write into the file
f= open("one.txt","w+")
#fi1 = open("stats_CB_120.txt","a") 
f.write('These are the statistics for t*wpi = 120 \n'
        #'H_mag_mean = %5.5f' %H_mag_mean + '\t  H_mag_std = %5.5f' %H_mag_std + '\t Hmag_threshold = %5.5f' %Hmag_threshold +'\n'  
        #'E_mag_mean = %5.5f' %E_mag_mean + '\t  E_mag_std = %5.5f' %E_mag_std + '\t Emag_threshold = %5.5f' %Emag_threshold +'\n'  
        #'J_mag_mean = %5.5f' %J_mag_mean + '\t  J_mag_std = %5.5f' %J_mag_std + '\t Jmag_threshold = %5.5f' %Jmag_threshold +'\n'  
        #'Vi_mag_mean = %5.5f' %Vi_mag_mean + '\t  Vi_mag_std = %5.5f' %Vi_mag_std + '\t Vi_threshold = %5.5f' %Vi_threshold +'\n'  
        #'Ve_mag_mean = %5.5f' %Ve_mag_mean + '\t  Ve_mag_std = %5.5f' %Ve_mag_std + '\t Ve_threshold = %5.5f' %Ve_threshold +'\n'  
        #'U_mag_mean = %5.5f' %U_mag_mean + '\t  U_mag_std = %5.5f' %U_mag_std + '\t U_threshold = %5.5f' %U_threshold +'\n' 
        'Ti_mean = %5.5f' %Ti_mean + '\t  Ti_std = %5.5f' %Ti_std + '\t Ti_threshold = %5.5f' %Ti_threshold +'\n'  
        'Te_mean = %5.5f' %Te_mean + '\t  Te_std = %5.5f' %Te_std + '\t Te_threshold = %5.5f' %Te_threshold +'\n'  
        'Vith_mean = %5.5f' %Vith_mean + '\t  Vith_std = %5.5f' %Vith_std + '\t Vith_threshold = %5.5f' %Vith_threshold +'\n'  
        'Veth_mean = %5.5f' %Veth_mean + '\t  Veth_std = %5.5f' %Veth_std + '\t Veth_threshold = %5.5f' %Veth_threshold +'\n'  
        'Ti_mean_R = %5.5f' %Ti_mean_R + '\t  Ti_std_R = %5.5f' %Ti_std_R + '\t Ti_threshold_R = %5.5f' %Ti_threshold_R +'\n'  
        'Te_mean_R = %5.5f' %Te_mean_R + '\t  Te_std_R = %5.5f' %Te_std_R + '\t Te_threshold_R = %5.5f' %Te_threshold_R +'\n'  
        'Vith_mean_R = %5.5f' %Vith_mean_R + '\t  Vith_std_R = %5.5f' %Vith_std_R + '\t Vith_threshold_R = %5.5f' %Vith_threshold_R +'\n'  
        'Veth_mean_R = %5.5f' %Veth_mean_R + '\t  Veth_std_R = %5.5f' %Veth_std_R + '\t Veth_threshold_R = %5.5f' %Veth_threshold_R +'\n'  
        #'E_par_mean = %5.5f' %E_par_mean + '\t  E_par_std = %5.5f' %E_par_std + '\t Epar_threshold = %5.5f' %Epar_threshold +'\n'  
        #'E_perp_mean = %5.5f' %E_perp_mean + '\t  E_perp_std = %5.5f' %E_perp_std + '\t Eperp_threshold = %5.5f' %Eperp_threshold +'\n'  
        #'J_par_mean = %5.5f' %J_par_mean + '\t  J_par_std = %5.5f' %J_par_std + '\t Jpar_threshold = %5.5f' %Jpar_threshold +'\n'  
        #'J_perp_mean = %5.5f' %J_perp_mean + '\t  J_perp_std = %5.5f' %J_perp_std + '\t Jperp_threshold = %5.5f' %Jperp_threshold +'\n'  
        )
f.close() 
# According to the conversation with george the places where the thermal energy increases
# are the palces where heating is taken place. Its important to check how the code is giving back
# the temperature tensor.  


# Calculate the histograms to see the distribution in order to identify acceleration sites
hist, binedges = np.histogram(J_mag, bins='auto', density=True)
hist_ve, binedges_ve = np.histogram(Ve_mag, bins='auto', density=True)
hist_vi, binedges_vi = np.histogram(Vi_mag, bins='auto', density=True)
hist_H, binedges_H = np.histogram(H_mag, bins='auto', density=True)
hist_Epar, binedges_Epar = np.histogram(E_par, bins='auto', density=True)


hist_Ti, binedges_Ti = np.histogram(Ti, bins='auto', density=True)
hist_Te, binedges_Te = np.histogram(Te, bins='auto', density=True)
hist_Vith, binedges_Vith = np.histogram(Vith, bins='auto', density=True)
hist_Veth, binedges_Veth = np.histogram(Veth, bins='auto', density=True)
hist_Ti_R, binedges_Ti_R = np.histogram(Ti_R, bins='auto', density=True)
hist_Te_R, binedges_Te_R = np.histogram(Te_R, bins='auto', density=True)
hist_Vith_R, binedges_Vith_R = np.histogram(Vith_R, bins='auto', density=True)
hist_Veth_R, binedges_Veth_R = np.histogram(Veth_R, bins='auto', density=True)


# Plot the histograms
f2, ax = plt.subplots()
#plt.plot(binedges[:-1],hist,label=r'$J_{Mag}$')
#plt.plot(binedges_ve[:-1],hist_ve,label=r'$Ve_{Mag}$')
#plt.plot(binedges_vi[:-1],hist_vi,label=r'$Vi_{Mag}$')
#plt.plot(binedges_H[:-1],hist_H,label=r'$H_{Mag}$')
#plt.plot(binedges_Epar[:-1],hist_Epar,label=r'$E_{par}$')
#plt.plot(binedges_vix[:-1],hist_vix,label=r'$vi_{x}$')
plt.plot(binedges_Ti[:-1],hist_Ti,label=r'$Ti$')
plt.plot(binedges_Te[:-1],hist_Te,label=r'$Te$')
#plt.plot(binedges_Vith[:-1],hist_Vith,label=r'$Vith$')
#plt.plot(binedges_Veth[:-1],hist_Veth,label=r'$Veth$')
plt.plot(binedges_Ti_R[:-1],hist_Ti_R,label=r'$Ti_R$')
plt.plot(binedges_Te_R[:-1],hist_Te_R,label=r'$Te_R$')
#plt.plot(binedges_Vith_R[:-1],hist_Vith_R,label=r'$Vith_R$')
#plt.plot(binedges_Veth_R[:-1],hist_Veth_R,label=r'$Veth_R$')
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.xlabel('Bins')
plt.ylabel('Probability')
#plt.title('Histogram normalised to the mean value')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(-0.05, 0.3)
#plt.ylim(0, 0.03)
#ax.axvline(x=96,color='k',linestyle='--')
plt.legend(loc='center right')
plt.grid(True)
plt.tight_layout()
plt.show()

f2, ax = plt.subplots()
ax.hist(r/np.mean(r), density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.show()


logvi=np.log10(Vi_mag[:,:,:])

hist_logvi, binedges_logvi = np.histogram(logvi,  bins='auto', density=True)


f2, ax = plt.subplots()
plt.plot(binedges_logvi[:-1],hist_logvi,label=r'$Vi_{Mag}$')
#plt.xlabel('Bins')
plt.ylabel('Probability')
#plt.title('Histogram of IQ')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(-0.05, 0.3)
#plt.ylim(0, 0.03)
#ax.axvline(x=96,color='k',linestyle='--')
plt.legend(loc='center right')
plt.grid(True)
plt.tight_layout()
plt.show()

f2.savefig("T_one_R_2.png", bbox_inches='tight')
###############################################################################################################


#------------------------------------------------------------------------------
#According to uritsky 2010
# Joule heating

Joule_1th_2_Uri=np.mean(magCurlB_2)+np.sqrt(np.mean(magCurlB_2**2)-(np.mean(magCurlB_2))**2)
J_sheet_ind_Uri = magCurlB_2 >  Joule_1th_2_Uri
J_sheet_2_Uri = magCurlB_2*J_sheet_ind_Uri

Joule_1th_2=np.mean(J_mag**2)+np.sqrt(np.mean(J_mag**4)-(np.mean(J_mag**2))**2) #It might be different using the current given by the simulation



J_mag_2=J_mag**2
J_sheet_ind = J_mag_2 >  Joule_1th_2
J_sheet_2 = J_mag_2*J_sheet_ind



# Kinetic dissipation
omega=magCurlvi**2
omega_1th_2=np.mean(magCurlvi**2)+np.sqrt(np.mean(magCurlvi**4)-(np.mean(magCurlvi**2))**2)
omega_sheet_ind_Uri = J_mag_2 >  Joule_1th_2_Uri
omega_sheet_2_Uri = J_mag_2*J_sheet_ind_Uri

J_mag_mean=np.mean(J_mag)


J_mag_std=np.std(J_mag)
J_sheet_ind = J_mag < J_mag_mean + J_mag_std 
#J_sheet_ind.astype(np.int) #false=0, true = 1
J_sheet = J_mag*J_sheet_ind

# Calculate the positions where the current is larger than a certain threshold









