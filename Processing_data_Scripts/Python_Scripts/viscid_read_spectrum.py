#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:09:46 2019
@author: jaa

This is the file that I will be using to make the plots henceforth
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


# This is to put the latex typography in the plots 
#--------------------------------------------------------------------
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# set the path to the file
project_dir ='/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/256_1_mode'
#project_dir ='/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256'
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

Utotal_av=[]; U_elec_av=[]; U_mag_av =[]; Ki_av =[]; Ke_av =[]; Ki_th_av =[]; Ke_th_av=[]; Current_av=[]; Hfield_av=[]
P3D1total=[]
for i in range(1, 50):    # This is a way to go throught the diferent times 
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
    U_elec = 0.5*Efield**2; U_mag = 0.5*Hfield**2
    Ki=0.5*vi**2; Ke=0.5*0.025*ve**2
    Ki_th=vi_th**2; Ke_th=ve_th**2 
    Utotal = U_elec + U_mag + Ki + Ke + Ki_th + Ke_th
    
    #average values over time
    num_cells= gx*gy*gz
    Utotal_av_i = np.sum(Utotal)/num_cells
    U_elec_av_i = np.sum(U_elec)/num_cells
    U_mag_av_i = np.sum(U_mag)/num_cells
    Ki_th_av_i=np.sum(Ki_th)/num_cells
    Ke_th_av_i=np.sum(Ke_th)/num_cells
    Ki_av_i=np.sum(Ki)/num_cells
    Ke_av_i=np.sum(Ke)/num_cells
    #Current_av_i=np.sum(Current)/num_cells
    
    Current_av_i=np.sqrt(np.mean(Current**2) - np.mean(Current)*np.mean(Current))
    Hfield_av_i=np.sqrt(np.mean(Hfield**2) - np.mean(Hfield)*np.mean(Hfield))
        
    Utotal_av.append(Utotal_av_i)
    U_elec_av.append(U_elec_av_i)
    U_mag_av.append(U_mag_av_i)
    Current_av.append(Current_av_i)
    Hfield_av.append(Hfield_av_i)
    Ki_th_av.append(Ki_th_av_i)
    Ke_th_av.append(Ke_th_av_i)
    Ki_av.append(Ki_av_i)
    Ke_av.append(Ke_av_i)

"""
    # Fourier transform
    #-------------------------------------------------------------------------
    # It must be calculated the fourier transform of each component of the vector field quantity 
    #f.activate_time(40) #Here I can change the step time 
    Fhx=np.fft.fftn(f['hx']); Fhy=np.fft.fftn(f['hy']); Fhz=np.fft.fftn(f['hz']);
    Fhx = np.fft.fftshift(Fhx); Fhy = np.fft.fftshift(Fhy); Fhz = np.fft.fftshift(Fhz) 
    CFhx=np.conjugate(Fhx); CFhy=np.conjugate(Fhy); CFhz=np.conjugate(Fhz)
    Bhat2=Fhx*CFhx + Fhy*CFhy + Fhz*CFhz     
    Bhat2_R_a=Bhat2.real

    # plt.pcolor(kxn,kyn,Bhat2.real[:,:,40])
    #This is the sum over row and columns
    P3D=[];    P3D1=[]
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
    
    P3D1total_inter=P3D1 # This is to create a large file with the spectrograms
    P3D1=[]
    P3D1total.append(P3D1total_inter)
    
P3D1total =  np.array(P3D1, dtype='float64')
P3D1total = np.transpose(P3D1)
"""

#-----------------------------------------------------------------------------
# Fourier transform for a single time step
#-------------------------------------------------------------------------
# It must be calculated the fourier transform of each component of the vector field quantity 
f.activate_time(1) #Here I can change the step time 
Fhx=np.fft.fftn(f['hx']); Fhy=np.fft.fftn(f['hy']); Fhz=np.fft.fftn(f['hz']);
Fhx = np.fft.fftshift(Fhx); Fhy = np.fft.fftshift(Fhy); Fhz = np.fft.fftshift(Fhz) 
CFhx=np.conjugate(Fhx); CFhy=np.conjugate(Fhy); CFhz=np.conjugate(Fhz)
Bhat2=Fhx*CFhx + Fhy*CFhy + Fhz*CFhz     

#--------------------------------------------------------------------------
#Calculate the k vectors which depends on the size of the system.
#Put the space resolution
di=1 
dx=0.078125; dy=0.078125; dz=0.078125
#grid point
gx=256; gy=256; gz=256
#Put the size of the box 
Lx=dx*gx*di; Ly=dy*gy*di; Lz=dz*gz*di
#for n in range(1,gx+1):
#kxn=n*2*np.pi/Lx
#kxn=(np.arange(gx))*(2*np.pi/Lx)    #This is what the fft return for the k axis

kxn=(np.arange(gx)-gx/2)*(2*np.pi/Lx)
kyn=(np.arange(gy)-gy/2)*(2*np.pi/Ly)
kzn=(np.arange(gz)-gz/2)*(2*np.pi/Lz) #257


Bhat2N=Bhat2/gx*gy*gz # Normilised by 
Bhat2_R_a=Bhat2.real


#This creates the meshgrid perpendicular to make the integration
#kxy= np.meshgrid(kxn,kyn)    
#Kper=[math.sqrt(x*x + y*y) for x, y in zip(kxn, kxn)]

kp=np.zeros((gx,gy))
for i in range(gx):
    for j in range(gy):
        kp[i,j]=math.sqrt(kxn[i]*kxn[i] + kyn[j]*kyn[j])

kp2=abs(kzn)
dk=2*np.pi/Lx
#dk=[]
#for i in range(1,len(kxn)):
#    dki=kxn[i]-kxn[i-1]
#    dk.append(dki)
Nkpar = int(np.max(kzn)/dk)
Nkper= int(np.max(kp)/dk)
#kpar=(0:nkpar-1).*dk;
kpar=(np.arange(Nkpar+1))*dk #+1 -1?
kper=(np.arange(Nkper+1))*dk

Ptemp=np.zeros((len(Bhat2[:,1,1]),Nkper)) #old kper vs. new kper
P2D=np.zeros((Nkpar,Nkper)) #(128,181)
NNkper=np.zeros((len(Bhat2[:,1,1]),Nkper))

i1=np.fix(kp/dk)
i2=i1+1
i3=np.fix(kp2/dk) #257
i1 =  np.array(i1, dtype='int')
i2 =  np.array(i2, dtype='int')
i3 =  np.array(i3, dtype='int')
#i3=i3[:-1]

k2=kp/dk - i1
k1=1-k2

for k in range(len(Bhat2[:,1,1])):
    P1=k1*Bhat2_R_a[:,:,k]
    P2=k2*Bhat2_R_a[:,:,k]
    for j in range(len(Bhat2[1,:,1])):
        for i in range(len(Bhat2[1,1,:])):
            if(i2[i,j] < Nkper):
                Ptemp[k,i2[i,j]]=Ptemp[k,i2[i,j]]+P2[i,j] # Interpolate power to next ring
                Ptemp[k,i1[i,j]]=Ptemp[k,i1[i,j]]+P1[i,j] # Interpolate power to previous ring
                #Update number of points interpolated for each ring:
                NNkper[i2[i,j]]=NNkper[i2[i,j]]+k2[i,j] 
                NNkper[i1[i,j]]=NNkper[i1[i,j]]+k1[i,j]

    if i3[k]<Nkpar:
        P2D[i3[k],:]=P2D[i3[k],:]+Ptemp[k,:]; #Fill array with kper values for each kpar

P2D=P2D/kper[2]; #why with kper[2]?

P2D=np.transpose(P2D)

kpar2=kpar[:-1]; 

#

P2D2=[]
for i in range(len(P2D[:,0])):
    aa=np.trim_zeros(P2D[:,1])
    P2D2.append(aa) #(127,126)
P2D2=np.array(P2D2,dtype='float64')


f2, ax = plt.subplots()
plt.contourf(kpar2,kper,np.log10(P2D), 20) #  kpar(126),   P2D(182,126)
#plt.contourf(kper,kpar2,np.log10(P2D2), 10)
ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title('2D Spectrum')
plt.ylabel(r"$K_{\perp} d_i$")
plt.xlabel(r"$K_{\|} d_i$")
plt.colorbar()
plt.show()



# To check the circles in the Kperp Kpar domain 
#--------------------------------------------------------------------------
"""
#circle1 = plt.Circle((4.5, 4.5), kp[0:4], color='b', fill=False)
circle = plt.Circle((0.63, 0.63), kp[0], color='b', fill=False)
circle1 = plt.Circle((0.63, 0.63), kp[1], color='b', fill=False)
circle2 = plt.Circle((0.63, 0.63), kp[2], color='b', fill=False)
circle3 = plt.Circle((0.63, 0.63), kp[3], color='b', fill=False)
circle4 = plt.Circle((0.63, 0.63), kp[4], color='b', fill=False)
index_2=np.where(kper < kp[1],kper, -1)
#index_3=np.where(index_2 > kp[0],index_2, -1)
fig, ax = plt.subplots() 
plt.pcolor(kxn,kyn,kper)
#plt.pcolor(kxn,kyn,index_2)
ax.add_artist(circle)
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(circle4)
plt.show()
#--------------------------------------------------------------------------"""

#Search the positions to integrate
#----------------------------------------------------------------------------
hola=[] 
for i in range(len(kp)): # len(kp)
    if i==0:
        newArr3 = np.where(np.less_equal(kper,kp[i]),kper, -1)
    else:
        newArr3 = np.where(np.logical_and(np.greater(kper,kp[i-1]),np.less_equal(kper,kp[i])), kper, -1)
    hola.append(newArr3)
#---------------------------------------------------------------

#------------------------------------------------------------

# plt.pcolor(kxn,kyn,Bhat2.real[:,:,40])
#This is the sum over row and columns
P3D=[];    P3D1=[]
for k in range(Bhat2_R_a.shape[2]):
    P1=k1[k]*Bhat2[:,:,k]
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

f2, ax = plt.subplots()
plt.contourf(kzn[129:],kp,np.log10(P3D1[:,128:]),20)
ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title('2D Spectrum')
plt.ylabel(r"$K_{\perp} d_i$")
plt.xlabel(r"$K_{\|} d_i$")
plt.colorbar()
plt.show()


plt.contourf(f['jz'][:,:,12],10)
plt.pcolor(f['jz'][:,:,12])

#f2, ax = plt.subplots()
#plt.contourf(np.log10(kzn[129:]),np.log10(kp),np.log10(P3D1[:,128:]),20)
#plt.title('2D Spectrum')
#plt.imshow(kzn,kp,np.log(np.abs(P3D1)))

"""
plt.pcolor(kxn,kyn,Bhat2.real[:,:,40])
plt.title(r'Magnetic field transform(k)')
plt.ylabel(r'$K_xn$')
plt.xlabel(r'$K_yn$')
plt.colorbar()
plt.show()
"""
#-----------------------------------------------------------------------------



# Plots of the time series
#-----------------------------------------------------------------------------
U_em_fields_av=[x + y for x, y in zip(U_elec_av, U_mag_av)] 


xticks=np.arange(49)*600 #This will change depending on the configuration

f3, ax = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.semilogy(xticks, Utotal_av, label=r'$Utotal_{av}$')
plt.semilogy(xticks, U_em_fields_av, label='U em fields')
plt.semilogy(xticks, Ki_th_av, label='$Ki_{th}$')
plt.semilogy(xticks, Ke_th_av, label='$Ke_{th}$')
plt.semilogy(xticks, Ki_av, label='$Ki$')
plt.semilogy(xticks, Ke_av, label='$Ke$')
plt.legend(loc='upper right')
plt.xlabel(r'$t \Omega_{pi}$')
ax.set_ylim([0.0001,5])
plt.show()

f4, ax = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.semilogy(xticks, Current_av, '-',label=r'$j_{rms}$')
plt.semilogy(xticks, Hfield_av, label=r'$B_{rms}$')
plt.legend(loc='upper right')
plt.xlabel(r'$t \Omega_{pi}$')
#ax.set_ylim([0.0001,1])
plt.show()




#-----------------------------------------------------------------------------
"""
#plt.pcolor(kzn,kp,P3D1)
f2, ax = plt.subplots()
plt.contourf(np.log10(kzn[129:]),np.log10(kp),np.log10(P3D1[:,128:]),20)
plt.pcolor(np.log10(kzn[128:]),np.log10(kp),np.log10(np.abs(P3D1[:,128:])))
#plt.imshow(kzn,kp,np.log(np.abs(P3D1)))
plt.title('2D Spectrum')
plt.ylabel('K_perp')
plt.xlabel('K_parallel')
plt.colorbar()
plt.show()

"""
