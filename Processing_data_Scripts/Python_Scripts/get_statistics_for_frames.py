#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:17:25 2020

@author: jaa
"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')


import numpy as np
import h5py
import glob, os

import gc
gc.collect()

#cd /run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256

#Call and read the file
###############################################################################


#Load the file
#for i in range(1,10):
#path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'
path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'
file=sorted(glob.glob(os.path.join(path, '*.h5')))
a=-1
#a=5
#Utotal_av=[]; U_elec_av=[]; U_mag_av =[]; Ki_av =[]; Ke_av =[]; Ki_th_av =[]; Ke_th_av=[]; J_rms_T=[]; B_rms_T=[]; Vi_rms_T=[]

#S_poynting=[]; E_elec_mag_V=[] 

for i in glob.glob( os.path.join(path, '*.h5')):
    a=a+1;   
#    print(file[a])   
    filename1= file[a]#'/disk/plasma2/jaa/CB8WAVES/CB8waves_04/pfd.002000_p000000.h5' 
    #Read the file
    print(filename1)
    h5_1st = h5py.File(filename1, 'r')
    #Get the fields
    ##############################################################################################################################################################
    e_field=h5_1st[list(h5_1st.keys())[16]] # Like this the outcome is a group
    h_field=h5_1st[list(h5_1st.keys())[17]] # ['hx', 'hy', 'hz']
    j_field=h5_1st[list(h5_1st.keys())[18]] # ['jx', 'jy', 'jz']
    v_field=h5_1st[list(h5_1st.keys())[24]] # ['vx_e', 'vx_i', 'vy_e', 'vy_i', 'vz_e', 'vz_i']
    #n_density=h5_1st[list(h5_1st.keys())[22]] #['n_e', 'n_i']
    #T_tensor=h5_1st[list(h5_1st.keys())[0]] #['Txx_e', 'Txx_i', 'Txy_e', 'Txy_i', 'Txz_e', 'Txz_i', 'Tyy_e', 'Tyy_i', 'Tyz_e', 'Tyz_i', 'Tzz_e', 'Tzz_i']
    
    # Electric field components
    #-----------------------------------------------------
    #Unreference the variables which are in a 3D
    
    ex=e_field['ex'] ; ex=ex['p0'] ; ex=ex['3d'] ; ex=np.array(ex) ; #ex=ex[:,:,:] ;# ex=np.transpose(ex, (1, 2, 0))
    ey=e_field['ey'] ; ey=ey['p0'] ; ey=ey['3d'] ; ey=np.array(ey) ;#ey=ey[:,:,:] ;# ey=np.transpose(ey, (1, 2, 0))
    ez=e_field['ez'] ; ez=ez['p0'] ; ez=ez['3d'] ; ez=np.array(ez) ;#ez=ez[:,:,:] ;# ez=np.transpose(ez, (1, 2, 0))
    
    # Magnetic field components
    #-----------------------------------------------------
    hx=h_field['hx'] ; hx=hx['p0'] ; hx=hx['3d'] ; hx=np.array(hx) ;#hx=hx[:,:,:] ;# hx=np.transpose(hx, (1, 2, 0))
    hy=h_field['hy'] ; hy=hy['p0'] ; hy=hy['3d'] ; hy=np.array(hy) ;#hy=hy[:,:,:] ;# hy=np.transpose(hy, (1, 2, 0))
    hz=h_field['hz'] ; hz=hz['p0'] ; hz=hz['3d'] ; hz=np.array(hz) ;#hz=hz[:,:,:] ;# hz=np.transpose(hz, (1, 2, 0))
    
    # Current field components
    #-----------------------------------------------------
    jx=j_field['jx'] ; jx=jx['p0'] ; jx=jx['3d'] ; jx=np.array(jx) ;#jx=jx[:,:,:] ;# jx=np.transpose(jx, (1, 2, 0))
    jy=j_field['jy'] ; jy=jy['p0'] ; jy=jy['3d'] ; jy=np.array(jy) ;#jy=jy[:,:,:] ;# jy=np.transpose(jy, (1, 2, 0))
    jz=j_field['jz'] ; jz=jz['p0'] ; jz=jz['3d'] ; jz=np.array(jz) ;#jz=jz[:,:,:] ;# jz=np.transpose(jz, (1, 2, 0))
    
    # Velocity Ions components
    #-----------------------------------------------------
    vx_i=v_field['vx_i'] ; vx_i=vx_i['p0'] ; vx_i=vx_i['3d'] ; vx_i=np.array(vx_i) ;#vx_i=vx_i[:,:,:] ;# vx_i=np.transpose(vx_i, (1, 2, 0))
    vy_i=v_field['vy_i'] ; vy_i=vy_i['p0'] ; vy_i=vy_i['3d'] ; vy_i=np.array(vy_i) ;#vy_i=vy_i[:,:,:] ;# vy_i=np.transpose(vy_i, (1, 2, 0))
    vz_i=v_field['vz_i'] ; vz_i=vz_i['p0'] ; vz_i=vz_i['3d'] ; vz_i=np.array(vz_i) ;#vz_i=vz_i[:,:,:] ;# vz_i=np.transpose(vz_i, (1, 2, 0))
    
    # Velocity electrons components
    #-----------------------------------------------------
    vx_e=v_field['vx_e'] ; vx_e=vx_e['p0'] ; vx_e=vx_e['3d'] ; vx_e=np.array(vx_e) ;#vx_e=vx_e[:,:,:] ;# vx_e=np.transpose(vx_e, (1, 2, 0))
    vy_e=v_field['vy_e'] ; vy_e=vy_e['p0'] ; vy_e=vy_e['3d'] ; vy_e=np.array(vy_e) ;#vy_e=vy_e[:,:,:] ;# vy_e=np.transpose(vy_e, (1, 2, 0))
    vz_e=v_field['vz_e'] ; vz_e=vz_e['p0'] ; vz_e=vz_e['3d'] ; vz_e=np.array(vz_e) ;#vz_e=vz_e[:,:,:] ;# vz_e=np.transpose(vz_e, (1, 2, 0))
    
    # Density Ions and electrons
    #-----------------------------------------------------
    #n_i=n_density['n_i'] ; n_i=n_i['p0'] ; n_i=n_i['3d'] ; n_i=n_i[:,:,:] ;# n_i=np.transpose(n_i, (1, 2, 0))
    #n_e=n_density['n_e'] ; n_e=n_e['p0'] ; n_e=n_e['3d'] ; n_e=n_e[:,:,:] ;# n_e=np.transpose(n_e, (1, 2, 0))
    
    # Temperature electron tensor components T_1st_single
    #-----------------------------------------------------
    #Txx_e = T_tensor['Txx_e'] ; Txx_e=Txx_e['p0'] ; Txx_e=Txx_e['3d'] ; Txx_e=Txx_e[:,:,:] ; #Txx_e=np.transpose(Txx_e, (1, 2, 0))
    #Txy_e = T_tensor['Txy_e'] ; Txy_e=Txy_e['p0'] ; Txy_e=Txy_e['3d'] ; Txy_e=Txy_e[:,:,:] ;#Txy_e=np.transpose(Txy_e, (1, 2, 0))
    #Txz_e = T_tensor['Txz_e'] ; Txz_e=Txz_e['p0'] ; Txz_e=Txz_e['3d'] ; Txz_e=Txz_e[:,:,:] ;#Txz_e=np.transpose(Txz_e, (1, 2, 0))
    #Tyy_e = T_tensor['Tyy_e'] ; Tyy_e=Tyy_e['p0'] ; Tyy_e=Tyy_e['3d'] ; Tyy_e=Tyy_e[:,:,:] ;#Tyy_e=np.transpose(Tyy_e, (1, 2, 0))
    #Tyz_e = T_tensor['Tyz_e'] ; Tyz_e=Tyz_e['p0'] ; Tyz_e=Tyz_e['3d'] ; Tyz_e=Tyz_e[:,:,:] ;#Tyz_e=np.transpose(Tyz_e, (1, 2, 0))
    #Tzz_e = T_tensor['Tzz_e'] ; Tzz_e=Tzz_e['p0'] ; Tzz_e=Tzz_e['3d'] ; Tzz_e=Tzz_e[:,:,:] ;#Tzz_e=np.transpose(Tzz_e, (1, 2, 0))
    # Temperature Ions tensor components
    #-----------------------------------------------------
    #Txx_i = T_tensor['Txx_i'] ; Txx_i=Txx_i['p0'] ; Txx_i=Txx_i['3d'] ; Txx_i=Txx_i[:,:,:] ;#Txx_i=np.transpose(Txx_i, (1, 2, 0))
    #Txy_i = T_tensor['Txy_i'] ; Txy_i=Txy_i['p0'] ; Txy_i=Txy_i['3d'] ; Txy_i=Txy_i[:,:,:] ;#Txy_i=np.transpose(Txy_i, (1, 2, 0))
    #Txz_i = T_tensor['Txz_i'] ; Txz_i=Txz_i['p0'] ; Txz_i=Txz_i['3d'] ; Txz_i=Txz_i[:,:,:] ;#Txz_i=np.transpose(Txz_i, (1, 2, 0))
    #Tyy_i = T_tensor['Tyy_i'] ; Tyy_i=Tyy_i['p0'] ; Tyy_i=Tyy_i['3d'] ; Tyy_i=Tyy_i[:,:,:] ;#Tyy_i=np.transpose(Tyy_i, (1, 2, 0))
    #Tyz_i = T_tensor['Tyz_i'] ; Tyz_i=Tyz_i['p0'] ; Tyz_i=Tyz_i['3d'] ; Tyz_i=Tyz_i[:,:,:] ;#Tyz_i=np.transpose(Tyz_i, (1, 2, 0))
    #Tzz_i = T_tensor['Tzz_i'] ; Tzz_i=Tzz_i['p0'] ; Tzz_i=Tzz_i['3d'] ; Tzz_i=Tzz_i[:,:,:] ;#Tzz_i=np.transpose(Tzz_i, (1, 2, 0))
    ##############################################################################################################################################################
    
    
    #Make operations with the loaded data
    ###############################################################################
    
    #plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension
    # ex.shape returns the dimensions of the matrix 
    """
    J_mag = np.sqrt(jx[:,:,:]**2 + jy[:,:,:]**2 + jz[:,:,:]**2)
    H_mag = np.sqrt(hx[:,:,:]**2 + hy[:,:,:]**2 + hz[:,:,:]**2)
    E_mag = np.sqrt(ex[:,:,:]**2 + ey[:,:,:]**2 + ez[:,:,:]**2)
    Vi_mag = np.sqrt(vx_i[:,:,:]**2 + vy_i[:,:,:]**2 + vz_i[:,:,:]**2)
    Ve_mag = np.sqrt(vx_e[:,:,:]**2 + vy_e[:,:,:]**2 + vz_e[:,:,:]**2)
    vi_th=np.sqrt((Txx_i[:,:,:] + Tyy_i[:,:,:] + Tzz_i[:,:,:])/3)
    ve_th=np.sqrt((Txx_e[:,:,:] + Tyy_e[:,:,:] + Tzz_e[:,:,:])/3)
    """
    mi_me=100;  
    J_mag = np.sqrt(jx**2 + jy**2 + jz**2)
    H_mag = np.sqrt(hx**2 + hy**2 + hz**2)
    E_mag = np.sqrt(ex**2 + ey**2 + ez**2)
    Vi_mag = np.sqrt(vx_i**2 + vy_i**2 + vz_i**2)
    Ve_mag = np.sqrt(vx_e**2 + vy_e**2 + vz_e**2)
    #Ti = (Txx_i + Tyy_i + Tzz_i)/3
    #Te = (Txx_e + Tyy_e + Tzz_e)/3
    #Vith = np.sqrt(2)*np.sqrt(np.abs(Ti))
    #Veth = np.sqrt(2*mi_me)*np.sqrt(np.abs(Te))   
    E_par = (ex*hx + ey*hy + ez*hz)/H_mag
    JE=ex*jx + ey*jy + ez*jz
    """
    H_mag_mean = np.mean(H_mag); H_mag_std = np.std(H_mag);
    E_mag_mean = np.mean(E_mag); E_mag_std = np.std(E_mag);
    J_mag_mean = np.mean(J_mag); J_mag_std = np.std(J_mag);
    Vi_mag_mean = np.mean(Vi_mag); Vi_mag_std = np.std(Vi_mag);
    Ve_mag_mean = np.mean(Ve_mag); Ve_mag_std = np.std(Ve_mag);
    Ti_mean = np.mean(Ti); Ti_std = np.std(Ti);
    Te_mean = np.mean(Te); Te_std = np.std(Te);
    Vith_mean = np.mean(Vith); Vith_std = np.std(Vith);
    Veth_mean = np.mean(Veth); Veth_std = np.std(Veth);
    E_par_mean = np.mean(E_par); E_par_std = np.std(E_par);
    JE_mean = np.mean(JE); JE_std = np.std(JE);
    """
    H_mag_rms = np.sqrt(np.mean(H_mag**2)); #H_mag_std = np.std(H_mag);
    E_mag_rms = np.sqrt(np.mean(E_mag**2)); #E_mag_std = np.std(E_mag);
    J_mag_rms = np.sqrt(np.mean(J_mag**2)); #J_mag_std = np.std(J_mag);
    Vi_mag_rms = np.sqrt(np.mean(Vi_mag**2)); #Vi_mag_std = np.std(Vi_mag);
    Ve_mag_rms = np.sqrt(np.mean(Ve_mag**2)); #Ve_mag_std = np.std(Ve_mag);
    #Ti_mean = np.mean(Ti); Ti_std = np.std(Ti);
    #Te_mean = np.mean(Te); Te_std = np.std(Te);
    #Vith_mean = np.mean(Vith); Vith_std = np.std(Vith);
    #Veth_mean = np.mean(Veth); Veth_std = np.std(Veth);
    E_par_rms = np.sqrt(np.mean(E_par**2)); #E_par_std = np.std(E_par);
    JE_rms = np.sqrt(np.mean(JE**2)); #JE_std = np.std(JE);
    
    
    #calculating thresholds    
    #H_mag_threshold = H_mag_mean + 4*H_mag_std
    #E_mag_threshold = E_mag_mean + 4*E_mag_std
    #J_mag_threshold = J_mag_mean + 4*J_mag_std
    #Vi_threshold = Vi_mag_mean + 4*Vi_mag_std
    #Ve_threshold = Ve_mag_mean + 4*Ve_mag_std
    #Ti_threshold = Ti_mean + 4*Ti_std
    #Te_threshold = Te_mean + 4*Te_std
    #Vith_threshold = Vith_mean + 4*Vith_std
    #Veth_threshold = Veth_mean + 4*Veth_std
    #E_par_threshold = E_par_mean + 4*E_par_std
    #JE_threshold = JE_mean + 4*JE_std
    
    
    #H_mag_thr1 = H_mag_mean + 3*H_mag_std
    #E_mag_thr1 = E_mag_mean + 3*E_mag_std
    #J_mag_thr1 = J_mag_mean + 3*J_mag_std
    #Vi_thr1 = Vi_mag_mean + 3*Vi_mag_std
    #Ve_thr1 = Ve_mag_mean + 3*Ve_mag_std
    #Ti_thr1 = Ti_mean + 3*Ti_std
    #Te_thr1 = Te_mean + 3*Te_std
    #Vith_thr1 = Vith_mean + 3*Vith_std
    #Veth_thr1 = Veth_mean + 3*Veth_std
    #E_par_thr1 = E_par_mean + 3*E_par_std
    #JE_thr1 = JE_mean + 3*JE_std
    
    ###################################################################################################
    os.chdir('/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics')  
    #os.chdir('/disk/plasmaz/jaa/Critical_Balance/Temperature_test/This_are_the_runs_2/one')
    # Write into the file
    
    name=str(a*24)
    f= open("%s_rms.txt" %name,"w+")
    f.write('These are the statistics for t*wpi = %s' %name +'\n'
        'H_mag_rms = %5.5f' %H_mag_rms + '\n'  
        'E_mag_rms = %5.5f' %E_mag_rms + '\n'  
        'J_mag_rms = %5.5f' %J_mag_rms + '\n'  
        'E_par_rms = %5.5f' %E_par_rms +'\n'
        'JE_rms = %5.5f' %JE_rms + '\n'
        'Vi_mag_rms = %5.5f' %Vi_mag_rms +'\n'  
        'Ve_mag_rms = %5.5f' %Ve_mag_rms +'\n'          
        #'H_mag_mean = %5.5f' %H_mag_mean + '\t  H_mag_std = %5.5f' %H_mag_std + '\t H_mag_thr1 = %5.5f' %H_mag_thr1 + '\t H_mag_threshold = %5.5f' %H_mag_threshold +'\n'  
        #'E_mag_mean = %5.5f' %E_mag_mean + '\t  E_mag_std = %5.5f' %E_mag_std + '\t E_mag_thr1 = %5.5f' %E_mag_thr1 + '\t E_mag_threshold = %5.5f' %E_mag_threshold +'\n'  
        #'J_mag_mean = %5.5f' %J_mag_mean + '\t  J_mag_std = %5.5f' %J_mag_std + '\t J_mag_thr1 = %5.5f' %J_mag_thr1 + '\t J_mag_threshold = %5.5f' %J_mag_threshold +'\n'  
        #'E_par_mean = %5.5f' %E_par_mean + '\t  E_par_std = %5.5f' %E_par_std + '\t E_par_thr1 = %5.5f' %E_par_thr1 + '\t E_par_threshold = %5.5f' %E_par_threshold +'\n'
        #'JE_mean = %5.5f' %JE_mean + '\t  JE_std = %5.5f' %JE_std + '\t JE_thr1 = %5.5f' %JE_thr1 + '\t JE_threshold = %5.5f' %JE_threshold +'\n'
        #'Vi_mag_mean = %5.5f' %Vi_mag_mean + '\t  Vi_mag_std = %5.5f' %Vi_mag_std + '\t Vi_thr1 = %5.5f' %Vi_thr1 + '\t Vi_threshold = %5.5f' %Vi_threshold +'\n'  
        #'Ve_mag_mean = %5.5f' %Ve_mag_mean + '\t  Ve_mag_std = %5.5f' %Ve_mag_std + '\t Ve_thr1 = %5.5f' %Ve_thr1 + '\t Ve_threshold = %5.5f' %Ve_threshold +'\n'  
     #   'Ti_mean = %5.5f' %Ti_mean + '\t  Ti_std = %5.5f' %Ti_std + '\t Ti_thr1 = %5.5f' %Ti_thr1 + '\t Ti_threshold = %5.5f' %Ti_threshold +'\n'  
     #   'Te_mean = %5.5f' %Te_mean + '\t  Te_std = %5.5f' %Te_std + '\t Te_thr1 = %5.5f' %Te_thr1 + '\t Te_threshold = %5.5f' %Te_threshold +'\n'  
     #   'Vith_mean = %5.5f' %Vith_mean + '\t  Vith_std = %5.5f' %Vith_std + '\t Vith_thr1 = %5.5f' %Vith_thr1 + '\t Vith_threshold = %5.5f' %Vith_threshold +'\n'  
     #   'Veth_mean = %5.5f' %Veth_mean + '\t  Veth_std = %5.5f' %Veth_std + '\t Veth_thr1 = %5.5f' %Veth_thr1 + '\t Veth_threshold = %5.5f' %Veth_threshold +'\n'  
        )
    f.close() 
    
    del ex, ey, ez, hx, hy, hz, jx, jy, jz, vx_i, vy_i, vz_i, vx_e, vy_e, vz_e 
    #del Txx_e, Tyy_e, Tzz_e, Txx_i, Tyy_i, Tzz_i, Ti, Te
    del J_mag, H_mag, E_mag, Vi_mag, Ve_mag, E_par, JE 
    #del  Vith, Veth




#############################################################################################################################
"""
# To read the file
#####################################################################
import numpy as np
import h5py
import glob, os
import gc
gc.collect()

####################################################################
f = open("0.txt", "r")
for x in f:
  lin = x
  y=lin.split('=')
  y=str(y[1])
  y=y.split('\t')
  print(y[0])
f.close()
#####################################################################

path = '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics'
file=sorted(glob.glob(os.path.join(path, '*.txt')))
a=-1

for i in glob.glob( os.path.join(path, '*.h5')):
    a=a+1;   
#    print(file[a])   
    filename2= file[a]#
    f = open(filename2, "r")
    lin = f.readline()
    lin = f.readline() #H_mag
    lin = f.readline() #E_mag 
    lin = f.readline() # J_mag
    lin = f.readline() # E_par
    lin = f.readline() # JE
    
    y=lin.split('=')
    mean=str(y[:][1]); std=str(y[:][2]); thr1=str(y[:][3]); threshold = y=str(y[:][4])
    mean=mean.split('\t')[0] ; std=std.split('\t')[0]; thr1=thr1.split('\t')[0]
    f.close()
###############################################################################
"""
#############################################################################################################################
    