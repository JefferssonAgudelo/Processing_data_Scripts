#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:09:46 2019

@author: jaa
"""


#Import modules
#---------------------------------------------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import viscid
from viscid.plot import vpyplot as vlt
import h5py
import pandas as pd
import math
viscid.calculator.evaluator.enabled = True


from scipy.fftpack import fft, ifft

from scipy import *
from scipy import integrate
from scipy.integrate import dblquad

#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# set the path to the file
#project_dir = '/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/1024_1st'
project_dir = '/disk/plasma4/jaa/PSC_things/trillian/TR'
#viscid.__name__
#viscid.__doc__
#viscid.__dict__
#The __dict__ attribute will return a dictionary object of module attributes, functions and other definitions and their respective values.
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Load the files fro all times
f = viscid.load_file(project_dir + "/pfd.xdmf", force_reload=True)
#With f.__dict__ it is possible to see which are the elements of the object
#These are the keys, all of them are 3D objects
# ['hx', 'hy', 'hz']
# ['jx', 'jy', 'jz']
# ['vx_e', 'vx_i', 'vy_e', 'vy_i', 'vz_e', 'vz_i']
#['n_e', 'n_i']
#['Txx_e', 'Txx_i', 'Txy_e', 'Txy_i', 'Txz_e', 'Txz_i', 'Tyy_e', 'Tyy_i', 'Tyz_e', 'Tyz_i', 'Tzz_e', 'Tzz_i']
#---------------------------------------------------------------------------


#Here make the operations
#---------------------------------------------------------------------------
#Defining the variables

Utotal_av=[]; U_elec_av=[]; U_mag_av =[]; Ki_av =[]; Ke_av =[]; Ki_th_av =[]; Ke_th_av=[]
for i in range(1, 80):    # This is a way to go throught the diferent times 
    f.activate_time(i)     
    
    # Calculating the magnitudes    
    Current=np.sqrt(f['jx']**2 + f['jy']**2 + f['jz']**2)
    Hfield=np.sqrt(f['hx']**2 + f['hy']**2 + f['hz']**2)
    Efield=np.sqrt(f['ex']**2 + f['ey']**2 + f['ez']**2)
    vi=np.sqrt(f['vx_i']**2 + f['vy_i']**2 + f['vz_i']**2)
    ve=np.sqrt(f['vx_e']**2 + f['vy_e']**2 + f['vz_e']**2)    
    vi_th=np.sqrt(f['Txx_i'] + f['Tyy_i'] + f['Tzz_i'])
    ve_th=np.sqrt(f['Txx_e'] + f['Tyy_e'] + f['Tzz_e'])

    #Energies here I have to multiply for the correct constants (FIXME) 
    U_elec = Efield**2; U_mag = Hfield**2
    Ki=vi**2; Ke=ve**2
    Ki_th=vi_th**2; Ke_th=ve_th**2 
    Utotal = U_elec + U_mag + Ki + Ke + Ki_th + Ke_th
    
    #average values over time
    num_cells= 16384
    Utotal_av_i = np.sum(Utotal)/num_cells
    U_elec_av_i = np.sum(U_elec)/num_cells
    U_mag_av_i = np.sum(U_mag)/num_cells
    
    Utotal_av.append(Utotal_av_i)
    U_elec_av.append(U_elec_av_i)
    U_mag_av.append(U_mag_av_i)


#It it easier to make them arrays
U_em_fields_av=[x + y for x, y in zip(U_elec_av, U_mag_av)] # sum the two lists
#U_em_fields_av_array=np.array(U_elec_av, dtype='float64') + np.array(U_mag_av, dtype='float64')


plt.semilogy(Utotal_av, label='Utotal_av')
plt.semilogy(U_em_fields_av, label='U_em_fields')
plt.legend(loc='upper right')
plt.show()


'''
The zip function is useful here, used with a list comprehension.
[x + y for x, y in zip(first, second)]
If you have a list of lists (instead of just two lists):
lists_of_lists = [[1, 2, 3], [4, 5, 6]]
[sum(x) for x in zip(*lists_of_lists)]
# -> [5, 7, 9]
'''
    
f.activate_time(40) 
Current=np.sqrt(f['jx']**2 + f['jy']**2 + f['jz']**2)
Hfield=np.sqrt(f['hx']**2 + f['hy']**2 + f['hz']**2)
Efield=np.sqrt(f['ex']**2 + f['ey']**2 + f['ez']**2)
vi=np.sqrt(f['vx_i']**2 + f['vy_i']**2 + f['vz_i']**2)
ve=np.sqrt(f['vx_e']**2 + f['vy_e']**2 + f['vz_e']**2)
vi_th=np.sqrt(f['Txx_i'] + f['Tyy_i'] + f['Tzz_i'])
ve_th=np.sqrt(f['Txx_e'] + f['Tyy_e'] + f['Tzz_e'])


U_ele=Efield**2; U_mag=Hfield**2
Ki=vi**2; Ke=ve**2
Ki_th=vi_th**2; Ke_th=ve_th**2 
    
Utotal = U_ele + U_mag + Ki + Ke + Ki_th + Ke_th

plt.semilogy(U_ele[4,4,:] , label='U_E')
plt.semilogy(U_mag[4,4,:], label='U_H')
plt.semilogy(Ki[4,4,:], label='Ki')
plt.semilogy(Ke[4,4,:], label='Ke')
plt.semilogy(Ki_th[4,4,:], label='Ki_th')
plt.semilogy(Ke_th[4,4,:], label='Ke_th')
plt.semilogy(Utotal[4,4,:], label='Utotal')
plt.legend(loc='upper right')
plt.show()

#plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
plt.show()    
    
    
    #Get the rms of the current normalise and get the positions
    Current_rms = np.sqrt(np.mean(Current**2))
    norm_current=Current/(Current_rms)

    newArr1 = Current[norm_current > 1]
    index_position_2=np.where(norm_current > 1)

    newArr2 = Current[Current > 6*Current_rms]
    index_position_2=np.where(Current > 6*Current_rms)
    plt.plot(Current[203,7,:], cmap="inferno")


    #Make the plots inside the for loop
    vlt.plot(Current[:,1,50], cmap="inferno")
    #plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
    plt.show()

   
#--------------------------------------------------------------------------
#Calculate the k vectors which depends on the size of the system.

#Put the space resolution
#dx=0.078125
di=1

dx=0.625; dy=0.625; dz=0.390625
#grid point
gx=8; gy=8; gz=256
#Put the size of the box 
Lx=dx*gx*di; Ly=dy*gy*di; Lz=dz*gz*di
#for n in range(1,gx+1):
    #kxn=n*2*np.pi/Lx

#kxn=(np.arange(gx))*(2*np.pi/Lx)    #This is what the fft return for the k axis

kxn=(np.arange(gx+1)-gx/2)*(2*np.pi/Lx)
kyn=(np.arange(gy+1)-gy/2)*(2*np.pi/Ly)
kzn=(np.arange(gz+1)-gz/2)*(2*np.pi/Lz)
#This creates the meshgrid perpendicular to make the integration
#kxy= np.meshgrid(kxn,kyn)    
#Kper=[math.sqrt(x*x + y*y) for x, y in zip(kxn, kxn)]

kperp=np.zeros((gx,gy))
for i in range(gx):
    for j in range(gy):
        kperp[i,j]=math.sqrt(kxn[i]*kxn[i] + kyn[j]*kyn[j])

dk=2*np.pi/Lx
Nk = int(np.max(kperp)/dk)
kp=(np.arange(Nk)+1)*dk

#circle1 = plt.Circle((4.5, 4.5), kp[0:4], color='b', fill=False)
circle = plt.Circle((0.63, 0.63), kp[0], color='b', fill=False)
circle1 = plt.Circle((0.63, 0.63), kp[1], color='b', fill=False)
circle2 = plt.Circle((0.63, 0.63), kp[2], color='b', fill=False)
circle3 = plt.Circle((0.63, 0.63), kp[3], color='b', fill=False)
circle4 = plt.Circle((0.63, 0.63), kp[4], color='b', fill=False)
#index_2=np.where(kperp < kp[1],kperp, -1)
#index_3=np.where(index_2 > kp[0],index_2, -1)
fig, ax = plt.subplots() 
plt.pcolor(kxn,kyn,kperp)
plt.pcolor(kxn,kyn,index_2)
ax.add_artist(circle)
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(circle4)
plt.show()


#Search the positions to integrate
#----------------------------------------------------------------------------
hola=[]
for i in range(len(kp)):
    if i==0:
        newArr3 = np.where(np.less_equal(kperp,kp[i]),kperp, -1)
    else:
        newArr3 = np.where(np.logical_and(np.greater(kperp,kp[i-1]),np.less_equal(kperp,kp[i])), kperp, -1)
    hola.append(newArr3)

# Fourier transform
#-------------------------------------------------------------------------
# It must be calculated the fourier transform of each component of the vector field quantity 

f.activate_time(40) #Here I can change the step time 
Fhx=np.fft.fftn(f['hx']); 
Fhy=np.fft.fftn(f['hy']);     
Fhz=np.fft.fftn(f['hz']);

Fhx = np.fft.fftshift(Fhx)
Fhy = np.fft.fftshift(Fhy)
Fhz = np.fft.fftshift(Fhz) 

CFhx=np.conjugate(Fhx)
CFhy=np.conjugate(Fhy)    
CFhz=np.conjugate(Fhz)

Bhat2=Fhx*CFhx + Fhy*CFhy + Fhz*CFhz     

#Bhat2_R_a=Bhat2.real
Bhat2_R_a=Bhat2.real

# plt.pcolor(kxn,kyn,Bhat2.real[:,:,40])

#This is the sum over row and columns
P3D=[]
P3D1=[]
for k in range(Bhat2_R_a.shape[2]):
    for i in range(len(kp)):
        hola2=hola[i]*Bhat2_R_a[:,:,k]
        P3D0_i= np.sum(hola2[hola2>=0])/len(hola2[hola2>=0]) #This is the normalize sum for the kperp
        #P3D0_i= np.sum((np.where(hola[i]*Bhat2_R_a[:,:,k] >= 0)))
        P3D.append(P3D0_i)
    P3D_inter=P3D
    P3D=[]
    P3D1.append(P3D_inter)
P3D1 =  np.array(P3D1, dtype='float64')
P3D1 = np.transpose(P3D1)

#plt.pcolor(kzn,kp,P3D1)
plt.pcolor(kzn,kp,np.log(np.abs(P3D1)))
#plt.imshow(kzn,kp,np.log(np.abs(P3D1)))
plt.title('2D Spectrum')
plt.ylabel('K_perp')
plt.xlabel('K_parallel')
plt.colorbar()
plt.show()



"""
#For non shifted fourier
#-------------------------------------------------------------------------
#For non shifted fourier
P3D=[]
P3D1=[]
for k in range(Bhat2_R_a.shape[2]):
    for i in range(Bhat2_R_a.shape[0]):
        P3D0_i= (np.sum(Bhat2_R_a[i,:,k][0:i+1]) + np.sum(Bhat2_R_a[:,i,k][0:i+1]) - Bhat2_R_a[i,i,k])/(2*i+1)
        P3D.append(P3D0_i)
    P3D_inter=P3D
    P3D=[]
    P3D1.append(P3D_inter)
P3D1 =  np.array(P3D1, dtype='float64')
P3D1 = np.transpose(P3D1)
plt.pcolor(P3D1)
plt.imshow(np.log(np.abs(P3D1)))
#plt.colorbar()
#plt.show()

hola4=[]
for k in range(256):
    for i in range(len(kp)):
        hola2=hola[i]*Bhat2[:,:,k]
        hola3 = np.sum((np.where(hola2.real >= 0)))
        hola3 = np.sum((np.where(hola[1]*Bhat2[:,:,1] >= 0)))
        hola4.append(hola3)
#index_2=np.where(kperp < kp[1],kperp, -1)
#index_3=np.where(index_2 > kp[0],index_2, -1)
#-------------------------------------------------------------------------
"""


#----------------------------
for i in range(Bhat2_R_a.shape[0]):
    kxn[i] = round(kxn[i], 2)
    kyn[i] = round(kyn[i], 2)
for i in range(Bhat2_R_a.shape[2]):
    kzn[i] = round(kzn[i], 2)

fig, ax = plt.subplots()
plt.pcolor(Fhx.real[:,:,40])
plt.xticks(np.arange(Bhat2_R_a.shape[0]),kxn, rotation=90)
plt.yticks(np.arange(8),kyn)
#plt.xticks(np.arange(256),kzn, rotation=90)
plt.show()


fshift = np.fft.fftshift(Fhx)
magnitude_spectrum = 20*np.log(np.abs(fshift))

plt.subplot(121),plt.imshow(f['hx'][:,3,:], cmap = 'inferno')
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(magnitude_spectrum[:,0,:], cmap = 'inferno')
plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()

plt.subplot(121),plt.imshow(magnitude_spectrum[:,1,:], cmap = 'inferno')
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(magnitude_spectrum[:,2,:], cmap = 'inferno')
plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()

plt.subplot(121),plt.imshow(magnitude_spectrum[:,3,:], cmap = 'inferno')
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(magnitude_spectrum[:,4,:], cmap = 'inferno')
plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()

#-----------------------------------------------------------------------------
freq0 = np.fft.fftfreq(f['hx'].shape[0])
freq1 = np.fft.fftfreq(f['hx'].shape[1])
freq2 = np.fft.fftfreq(f['hx'].shape[2])
l=np.arange(1,9)
ind=[math.sqrt(2*x*x) for x, y in zip(l, l)]

data = np.ones((4, 4, 4))
data_array=np.array(data, dtype='float64')
dbtt=[]
dbt2=[]
dbt3=[]
for k in range(4):        
    for i in range(4):
        dbt= np.sum(data_array[i,:,k][0:i+1]) + np.sum(data_array[:,i,k][0:i+1]) - data_array[i,i,k] 
        dbtt.append(dbt)    
    dbt2=dbtt
    dbtt=[]
    dbt3.append(dbt2)   
dbt3 =  np.array(dbt3, dtype='float64')    
   

db3=data_array[2].sum()
db3c=data_array[:,2].sum()
dbt= data_array[i].sum()+ data_array[:,i].sum() - data_array[i,i]
np.sum(data_array[2][0:nn-1]) + np.sum(data_array[:,2][0:nn-1]) - data_array[nn-2,nn-2]


plt.plot(freq0, Bhat2[:,4,40].real, '*')
plt.pcolor(Bhat2.real[:,:,40])
plt.imshow(np.log(np.abs(Bhat2[:,:,40])))

bfield=np.sqrt(f['hx']**2 + f['hy']**2 + f['hz']**2)
plt.imshow(np.log(np.abs(bfield[4,:,:])))

#plt.pcolor(np.sum(bfield,axis=0)) #integration over the first dimension
#plt.colorbar()
#plt.show()
#-----------------------------------------------------------------------------