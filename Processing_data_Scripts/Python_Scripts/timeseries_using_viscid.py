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

import csv
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# set the path to the file
#project_dir = '/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/1024_1st'
#project_dir = '/disk/plasma4/jaa/PSC_things/trillian/TR'
project_dir = '/scratch/dp126/dc-agud1/JOBS_SCRATCH/TR_RUNS/TR_3/WEAK_SCALING/1024_1st'
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
for i in range(1, 44):    # This is a way to go throught the diferent times 
    f.activate_time(i)     
    
    # Calculating the magnitudes    
    Current=np.sqrt(f['jx']**2 + f['jy']**2 + f['jz']**2)
    Hfield=np.sqrt(f['hx']**2 + f['hy']**2 + f['hz']**2)
    Efield=np.sqrt(f['ex']**2 + f['ey']**2 + f['ez']**2)
    vi=np.sqrt(f['vx_i']**2 + f['vy_i']**2 + f['vz_i']**2)
    ve=np.sqrt(f['vx_e']**2 + f['vy_e']**2 + f['vz_e']**2)    
    vi_th=np.sqrt(f['Txx_i'] + f['Tyy_i'] + f['Tzz_i'])
    ve_th=np.sqrt(f['Txx_e'] + f['Tyy_e'] + f['Tzz_e'])

    #Energies
    U_elec = Efield**2; U_mag = Hfield**2
    Ki=vi**2; Ke=ve**2
    Ki_th=vi_th**2; Ke_th=ve_th**2 
    Utotal = U_elec + U_mag + Ki + Ke + Ki_th + Ke_th
    
    #average values over time
    #num_cells= 16384
    num_cells= 67108864
    Utotal_av_i = np.sum(Utotal)/num_cells
    U_elec_av_i = np.sum(U_elec)/num_cells
    U_mag_av_i = np.sum(U_mag)/num_cells
    
    Utotal_av.append(Utotal_av_i)
    U_elec_av.append(U_elec_av_i)
    U_mag_av.append(U_mag_av_i)

U_em_fields_av=[x + y for x, y in zip(U_elec_av, U_mag_av)] # sum the two lists


combined = np.vstack((U_em_fields_av, Utotal_av)).T
np.savetxt("time_series_averages.csv", combined , delimiter=",", header=' U_em, Ut  ')




#--------------------------------------------------------------------------------
# Open and read the file that I bring from the cluster
with open('/disk/plasma4/jaa/PSC_things/Paraview_Jeff/Python_things/time_series_averages.csv', 'r') as csv_file_line:    
    csv_reader_line = csv.reader(csv_file_line, delimiter=',')

    list_hx_b=[]
    list_hy_b=[]
    line_count = 0
    
    for row in csv_reader_line:
            list_hx_b.append(row[0])
            list_hy_b.append(row[1])
    
    del list_hx_b[0]
    del list_hy_b[0]

array_U_em=np.array(list_hx_b, dtype='float64')
array_U_t=np.array(list_hy_b, dtype='float64')


l=np.arange(43)
ll=600*l

f2, ax = plt.subplots()
plt.semilogy(ll, array_U_t,label='Total energy density')
plt.semilogy(ll, array_U_em,label='EM energy density')
plt.legend(loc='upper right')
plt.title('Non dimensional Average Energy density time series')
plt.ylabel('Energy density')
plt.xlabel('Time steps')
plt.tight_layout()
plt.show()
f2.savefig("Energy_density_8_waves.png", bbox_inches='tight')

#------------------------------------------------------------------------------