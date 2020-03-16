#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:06:02 2019

@author: jaa
"""

import numpy as np
#import math
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
#import csv
#from scipy import signal
#from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import h5py
import glob, os
#cd /run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256

#Call and read the file
###############################################################################


#Load the file
#for i in range(1,10):
#path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'
path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1'

file=sorted(glob.glob(os.path.join(path, '*.h5')))
a=-1
#Utotal_av=[]; U_elec_av=[]; U_mag_av =[]; Ki_av =[]; Ke_av =[]; Ki_th_av =[]; Ke_th_av=[]; J_rms_T=[]; B_rms_T=[]; Vi_rms_T=[]

S_poynting=[]; E_elec_mag_V=[] 

for i in glob.glob( os.path.join(path, '*.h5')):
    a=a+1;   
#    print(file[a])   
    filename1= file[a]#'/disk/plasma2/jaa/CB8WAVES/CB8waves_04/pfd.002000_p000000.h5' 
    #Read the file
    print(filename1)
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
    
    ex=e_field['ex'] ; ex=ex['p0'] ; ex=ex['3d'] ; ex=ex[:,:,:] ; ex=np.transpose(ex, (1, 2, 0))
    ey=e_field['ey'] ; ey=ey['p0'] ; ey=ey['3d'] ; ey=ey[:,:,:] ; ey=np.transpose(ey, (1, 2, 0))
    ez=e_field['ez'] ; ez=ez['p0'] ; ez=ez['3d'] ; ez=ez[:,:,:] ; ez=np.transpose(ez, (1, 2, 0))
    
    # Magnetic field components
    #-----------------------------------------------------
    hx=h_field['hx'] ; hx=hx['p0'] ; hx=hx['3d'] ; hx=hx[:,:,:] ; hx=np.transpose(hx, (1, 2, 0))
    hy=h_field['hy'] ; hy=hy['p0'] ; hy=hy['3d'] ; hy=hy[:,:,:] ; hy=np.transpose(hy, (1, 2, 0))
    hz=h_field['hz'] ; hz=hz['p0'] ; hz=hz['3d'] ; hz=hz[:,:,:] ; hz=np.transpose(hz, (1, 2, 0))
    
    # Current field components
    #-----------------------------------------------------
    jx=j_field['jx'] ; jx=jx['p0'] ; jx=jx['3d'] ; jx=jx[:,:,:] ; jx=np.transpose(jx, (1, 2, 0))
    jy=j_field['jy'] ; jy=jy['p0'] ; jy=jy['3d'] ; jy=jy[:,:,:] ; jy=np.transpose(jy, (1, 2, 0))
    jz=j_field['jz'] ; jz=jz['p0'] ; jz=jz['3d'] ; jz=jz[:,:,:] ; jz=np.transpose(jz, (1, 2, 0))
    
    # Velocity Ions components
    #-----------------------------------------------------
    #vx_i=v_field['vx_i'] ; vx_i=vx_i['p0'] ; vx_i=vx_i['3d'] #; vx_i=vx_i[:,:,:] ; vx_i=np.transpose(vx_i, (1, 2, 0))
    #vy_i=v_field['vy_i'] ; vy_i=vy_i['p0'] ; vy_i=vy_i['3d'] #; vy_i=vy_i[:,:,:] ; vy_i=np.transpose(vy_i, (1, 2, 0))
    #vz_i=v_field['vz_i'] ; vz_i=vz_i['p0'] ; vz_i=vz_i['3d'] #; vz_i=vz_i[:,:,:] ; vz_i=np.transpose(vz_i, (1, 2, 0))
    
    # Velocity electrons components
    #-----------------------------------------------------
    #vx_e=v_field['vx_e'] ; vx_e=vx_e['p0'] ; vx_e=vx_e['3d'] #; vx_e=vx_e[:,:,:] ; vx_e=np.transpose(vx_e, (1, 2, 0))
    #vy_e=v_field['vy_e'] ; vy_e=vy_e['p0'] ; vy_e=vy_e['3d'] #; vy_e=vy_e[:,:,:] ; vy_e=np.transpose(vy_e, (1, 2, 0))
    #vz_e=v_field['vz_e'] ; vz_e=vz_e['p0'] ; vz_e=vz_e['3d'] #; vz_e=vz_e[:,:,:] ; vz_e=np.transpose(vz_e, (1, 2, 0))
    
    # Density Ions and electrons
    #-----------------------------------------------------
    #n_i=n_density['n_i'] ; n_i=n_i['p0'] ; n_i=n_i['3d'] ; #n_i=n_i[:,:,:] ; n_i=np.transpose(n_i, (1, 2, 0))
    #n_e=n_density['n_e'] ; n_e=n_e['p0'] ; n_e=n_e['3d'] ; #n_e=n_e[:,:,:] ; n_e=np.transpose(n_e, (1, 2, 0))
    
    
    # Temperature electron tensor components T_1st_single
    #-----------------------------------------------------
    #Txx_e = T_tensor['Txx_e'] ; Txx_e=Txx_e['p0'] ; Txx_e=Txx_e['3d'] ; #Txx_e=np.transpose(Txx_e, (1, 2, 0))
    #Txy_e = T_tensor['Txy_e'] ; Txy_e=Txy_e['p0'] ; Txy_e=Txy_e['3d'] ; #Txy_e=np.transpose(Txy_e, (1, 2, 0))
    #Txz_e = T_tensor['Txz_e'] ; Txz_e=Txz_e['p0'] ; Txz_e=Txz_e['3d'] ; #Txz_e=np.transpose(Txz_e, (1, 2, 0))
    #Tyy_e = T_tensor['Tyy_e'] ; Tyy_e=Tyy_e['p0'] ; Tyy_e=Tyy_e['3d'] ; #Tyy_e=np.transpose(Tyy_e, (1, 2, 0))
    #Tyz_e = T_tensor['Tyz_e'] ; Tyz_e=Tyz_e['p0'] ; Tyz_e=Tyz_e['3d'] ; #Tyz_e=np.transpose(Tyz_e, (1, 2, 0))
    #Tzz_e = T_tensor['Tzz_e'] ; Tzz_e=Tzz_e['p0'] ; Tzz_e=Tzz_e['3d'] ; #Tzz_e=np.transpose(Tzz_e, (1, 2, 0))
    # Temperature Ions tensor components
    #-----------------------------------------------------
    #Txx_i = T_tensor['Txx_i'] ; Txx_i=Txx_i['p0'] ; Txx_i=Txx_i['3d'] ; #Txx_i=np.transpose(Txx_i, (1, 2, 0))
    #Txy_i = T_tensor['Txy_i'] ; Txy_i=Txy_i['p0'] ; Txy_i=Txy_i['3d'] ; #Txy_i=np.transpose(Txy_i, (1, 2, 0))
    #Txz_i = T_tensor['Txz_i'] ; Txz_i=Txz_i['p0'] ; Txz_i=Txz_i['3d'] ; #Txz_i=np.transpose(Txz_i, (1, 2, 0))
    #Tyy_i = T_tensor['Tyy_i'] ; Tyy_i=Tyy_i['p0'] ; Tyy_i=Tyy_i['3d'] ; #Tyy_i=np.transpose(Tyy_i, (1, 2, 0))
    #Tyz_i = T_tensor['Tyz_i'] ; Tyz_i=Tyz_i['p0'] ; Tyz_i=Tyz_i['3d'] ; #Tyz_i=np.transpose(Tyz_i, (1, 2, 0))
    #Tzz_i = T_tensor['Tzz_i'] ; Tzz_i=Tzz_i['p0'] ; Tzz_i=Tzz_i['3d'] ; #Tzz_i=np.transpose(Tzz_i, (1, 2, 0))
    ###############################################################################
    
    
    
    """
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
    """

    #Make operations with the loaded data
    ###############################################################################
    
    #plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
    # ex.shape returns the dimensions of the matrix

    """
    mi_me=100;    
    J_mag = np.sqrt(jx[:,:,:]**2 + jy[:,:,:]**2 + jz[:,:,:]**2)
    H_mag = np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2 + hz[:,:,:]**2)
    E_mag = np.sqrt(ex[:,:,:]**2 + ey[:,:,:]**2 + ez[:,:,:]**2)
    Vi_mag = np.sqrt(vx_i[:,:,:]**2 + vy_i[:,:,:]**2 + vz_i[:,:,:]**2)
    Ve_mag = np.sqrt(vx_e[:,:,:]**2 + vy_e[:,:,:]**2 + vz_e[:,:,:]**2)
    vi_th=np.sqrt((Txx_i[:,:,:] + Tyy_i[:,:,:] + Tzz_i[:,:,:])/3)
    ve_th=np.sqrt((Txx_e[:,:,:] + Tyy_e[:,:,:] + Tzz_e[:,:,:])/3)
    """
    
    # Calculating the Poyting theorem conservation
    nx=len(hx[:,0,0])-1 # This takes a lot of time better just put by hand
    ny=len(hx[0,:,0])-1
    nz=len(hx[0,0,:])-1

    # The magnetic energy density integrated over the whole box
    H_mag = np.sqrt(hx*hx + hy*hy* + hz*hz)
    E_mag = np.sqrt(ex*ex + ey*ey* + ez*ez)
    E_elec_mag_V_i = np.sum(0.5*(H_mag**2 + E_mag**2))    
    
    # Calculate the integral of the over the boundary
    ExBda_x = np.sum(ey[nx,:,:]*hz[nx,:,:] - ez[nx,:,:]*hy[nx,:,:]) - np.sum(ey[0,:,:]*hz[0,:,:] - ez[0,:,:]*hy[0,:,:]) 
    ExBda_y = np.sum(ez[:,ny,:]*hx[:,ny,:] - ex[:,ny,:]*hz[:,ny,:]) - np.sum(ez[:,0,:]*hx[:,0,:] - ex[:,0,:]*hz[:,0,:]) 
    ExBda_z = np.sum(ex[:,:,nz]*hy[:,:,nz] - ey[:,:,nz]*hx[:,:,nz]) - np.sum(ex[:,:,0]*hy[:,:,0] - ey[:,:,0]*hx[:,:,0]) 

    ExB_da = ExBda_x + ExBda_y + ExBda_z
    EdotJ_V = np.sum(ex*jx + ey*jy + ez*jz)
    S_poynting_i = ExB_da + EdotJ_V

    """
    U_mag =  (n_e[:,:,:]*Ve_mag[:,:,:] + n_i[:,:,:]*Vi_mag[:,:,:]*mi_me)/(n_e[:,:,:] + n_i[:,:,:]*mi_me)
    num_cells=2016*400*400
    Jrms_i=np.sqrt((np.sum(J_mag[:,:,:]**2)/num_cells)-(np.sum(J_mag[:,:,:])/num_cells)**2)
    Brms_i=np.sqrt((np.sum(H_mag[:,:,:]**2)/num_cells)-(np.sum(H_mag[:,:,:])/num_cells)**2)
    Virms_i=np.sqrt((np.sum(Vi_mag[:,:,:]**2)/num_cells)-(np.sum(Vi_mag[:,:,:])/num_cells)**2)    
    
    #clear variables
        U_elec = 0.5*np.array(E_mag[:,:,:]**2); U_mag = 0.5*np.array(H_mag[:,:,:]**2)
    Ki=0.5*np.array(Vi_mag[:,:,:]**2); Ke=(0.5/mi_me)*np.array(Ve_mag[:,:,:]**2)
    Ki_th=0.5*np.array(vi_th[:,:,:]**2); Ke_th=(0.5/mi_me)*np.array(ve_th[:,:,:]**2) 
    Utotal = U_elec[:,:,:] + U_mag[:,:,:] + Ki[:,:,:] + Ke[:,:,:] + Ki_th[:,:,:] + Ke_th[:,:,:]
    
    Utotal_av_i = np.sum(Utotal[:,:,:])/num_cells
    U_elec_av_i = np.sum(U_elec[:,:,:])/num_cells
    U_mag_av_i = np.sum(U_mag[:,:,:])/num_cells    
    
    J_rms_T.append(Jrms_i) 
    B_rms_T.append(Brms_i) 
    Vi_rms_T.append(Virms_i)
    Utotal_av.append(Utotal_av_i)
    U_elec_av.append(U_elec_av_i)
    U_mag_av.append(U_mag_av_i)
    """

    S_poynting.append(S_poynting_i)    
    E_elec_mag_V.append(E_elec_mag_V_i)

S_poynting = np.array(S_poynting)
E_elec_mag_V = np.array(E_elec_mag_V)

dEem_dt=[]
for i in range(len(E_elec_mag_V)-1):
    dEem_dt_i = E_elec_mag_V[i+1] - E_elec_mag_V[i]
    dEem_dt.append(dEem_dt_i)

PConserv = dEem_dt + S_poynting[:-1]

#U_em_fields_av=[x + y for x, y in zip(U_elec_av, U_mag_av)]


#Make the plot
#####################################################################
#ll=len(J_rms_T)
ll=len(PConserv)
l=np.arange(ll)
ll=24*l #400*0.06


"""
f2, ax = plt.subplots()
plt.plot(ll, J_rms_T,label=r'$J^{rms}$')
plt.plot(ll, B_rms_T,label=r'$B^{rms}$')
plt.plot(ll, Vi_rms_T,label=r'$V_{i}^{rms}$')
ax.axvline(x=96,color='k',linestyle='--')
plt.legend(loc='center right')
#plt.title('Non dimensional Average Energy density time series')
#plt.ylabel('Energy density')
plt.xlabel(r'$\omega_{pi} t$')
plt.tight_layout()
plt.show()
f3, ax = plt.subplots()
plt.semilogy(ll, Utotal_av,label=r'$\langle U_{T} \rangle$')
plt.semilogy(ll, U_elec_av,label=r'$\langle U_{E} \rangle$')
plt.semilogy(ll, U_mag_av,label=r'$\langle U_{B} \rangle$')
plt.legend(loc='center right')
#plt.title('Non dimensional Average Energy density time series')
#plt.ylabel('Energy density')
plt.xlabel(r'$\omega_{pi} t$')
plt.tight_layout()
plt.show()
f4, ax = plt.subplots()
plt.plot(ll, Utotal_av,label=r'$\langle U_{T} \rangle$')
plt.plot(ll, U_em_fields_av,label=r'$\langle U_{EM} \rangle$')
plt.legend(loc='center right')
#plt.title('Non dimensional Average Energy density time series')
#plt.ylabel('Energy density')
plt.xlabel(r'$\omega_{pi} t$')
plt.tight_layout()
plt.show()
f2.savefig("JBVi_rms_2.png", bbox_inches='tight')
f3.savefig("Energy_density_8_waves_2.png", bbox_inches='tight')
f4.savefig("EM_UT_8_waves_2.png", bbox_inches='tight')
"""

f5, ax = plt.subplots()
plt.plot(ll, PConserv)
plt.legend(loc='center right')
plt.title('Energy Conservation')
#plt.ylabel('Energy Conservation')
plt.xlabel(r'$\omega_{pi} t$')
plt.tight_layout()
plt.show()
f5.savefig("poynting_time.png", bbox_inches='tight')


#cd /disk/plasma2/jaa/CB8WAVES/CB8waves_04/



"""
    E_par = (ex[:,:,:]*hx[:,:,:] + ey[:,:,:]*hy[:,:,:] + ez[:,:,:]*hz[:,:,:])/H_mag
    theta_e = np.arccos((ex[:,:,:]*hx[:,:,:] + ey[:,:,:]*hy[:,:,:] + ez[:,:,:]*hz[:,:,:])/(E_mag*H_mag))
    E_perp = E_mag*np.sin(theta_e)
    
    J_par = (jx[:,:,:]*hx[:,:,:] + jy[:,:,:]*hy[:,:,:] + jz[:,:,:]*hz[:,:,:])/H_mag
    theta_j = np.arccos((jx[:,:,:]*hx[:,:,:] + jy[:,:,:]*hy[:,:,:] + jz[:,:,:]*hz[:,:,:])/(J_mag*H_mag))
    J_perp = J_mag*np.sin(theta_j)
    
    JE_mag=ex[:,:,:]*jx[:,:,:] + ey[:,:,:]*jy[:,:,:] + ez[:,:,:]*jz[:,:,:]
"""









"""
#get the maximum
    J_mag_max=np.amax(J_mag[:,:,:])
    index_position=np.where(J_mag == J_mag_max)    
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
"""

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

