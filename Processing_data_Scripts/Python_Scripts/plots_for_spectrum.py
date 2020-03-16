#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:48:21 2019

@author: jaa
"""

import math
import numpy as np
import csv
from scipy import signal
from scipy.stats import linregress

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter  # useful for `logit` scale

import viscid
import pandas as pd
viscid.calculator.evaluator.enabled = True
from viscid.plot import vpyplot as vlt

from matplotlib.backends.backend_pdf import PdfPages


# This load the variables readed by 
import h5_filereader_jeff as reader
print(reader.ex.shape)

ex=reader.ex


#plt.pcolor(ex[65,:,:]) #plot in the 65 (middle)

#plt.pcolor(np.sum(ex,axis=0)) #integration over the first dimension (this dimension is the z direction)



# Read file from a overline file 
############################################################
with open('/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/1024_stuff/1st/1024_4800_1st_line.csv', 'r') as csv_file_line:    
    csv_reader_line = csv.reader(csv_file_line, delimiter=',')
    
    list_hx_b=[]
    list_hy_b=[]
    list_hz_b=[]
    line_count = 0
    
    for row in csv_reader_line:
            list_hx_b.append(row[3])
            list_hy_b.append(row[4])
            list_hz_b.append(row[5])    

    del list_hx_b[0]
    del list_hy_b[0]
    del list_hz_b[0]

array_hx_b=np.array(list_hx_b, dtype='float64')
array_hy_b=np.array(list_hy_b, dtype='float64')
array_hz_b=np.array(list_hz_b, dtype='float64')

array_BT_b=np.sqrt(np.power(array_hx_b,2)+np.power(array_hy_b,2)+np.power(array_hz_b,2))


# Generate a chirp signal
############################################################
sig_BT_b = array_hx_b#array_BT_b 
sig_BT2_b = np.sqrt(np.power((array_BT_b -array_BT_b.mean()),2)) #The variations




wave_n_b, lenght_b, spectrogram_b = signal.spectrogram(sig_BT_b)
wave_n_2b, lenght_2b, spectrogram_2b = signal.spectrogram(sig_BT2_b)


plt.figure(figsize=(5, 4))
plt.imshow(spectrogram_b, aspect='auto', cmap='hot_r', origin='lower')
plt.title('Magnetic field')
plt.ylabel('Wavenumber band')
plt.xlabel('Space window')
plt.tight_layout()
plt.show()

plt.figure(figsize=(5, 4))
plt.imshow(spectrogram_2b, aspect='auto', cmap='hot_r', origin='lower')
plt.title('Variations of Bfield')
plt.ylabel('Wavenumber band')
plt.xlabel('Space window')
plt.tight_layout()
plt.show()

wave_n_b, psd_b = signal.welch(sig_BT_b)
plt.figure(figsize=(5, 4))
plt.semilogx(wave_n_b, psd_b)
plt.title('PSD: power spectral density hx')
plt.yscale('log')
plt.xlabel('Wavenumber')
plt.ylabel('Power')
plt.tight_layout()
plt.show()


wave_n_2b, psd_2b = signal.welch(sig_BT2_b )
plt.figure(figsize=(5, 4))
plt.semilogx(wave_n_2b, psd_2b)
plt.title('PSD: power spectral density var')
plt.yscale('log')
plt.xlabel('Wavenumber')
plt.ylabel('Power')
plt.tight_layout()
plt.show()




#Regression
############################################################
log_wave_n_b=np.log(wave_n_b)

#To remove nan and inf values
ii_b=np.isfinite(log_wave_n_b)
log2_b=log_wave_n_b[ii_b]
psd2_b=psd_b[1:] # Remove the first value 

slope_b, intercept_b, r_value_b, p_value_b, std_err_b = linregress(log2_b, psd2_b)
print(slope_b)
psd3_b=slope_b*log2_b + intercept_b

"""
fig2=plt.figure(1)
plt.semilogx(wave_n_b, psd_b)
plt.semilogx(np.exp(log2_b), psd3_b, 'r')
plt.yscale('log')
plt.title('PSD: power spectral density')
#plt.text(60, .025, r'$\alpha_k='%slope)
plt.xlabel('Wavenumber ($kd_{i}$)')
plt.ylabel('Power')
plt.tight_layout()
plt.show(fig2)
"""
############################################################

dt = 0.0005
t= array_hx_b

x1 = sig_BT_b   # the signal
x2= sig_BT2_b
NFFT = 1024  # the length of the windowing segments
Fs = int(1.0 / dt)  # the sampling frequency

"""
plt.figure(figsize=(5, 4))
Pxx, freqs, bins, im = plt.specgram(x1, NFFT=NFFT, Fs=Fs, noverlap=900)
plt.show()
"""

fig, (ax1, ax2) = plt.subplots(nrows=2)
Pxx, freqs, bins, im = ax1.specgram(x1, NFFT=NFFT, Fs=Fs, noverlap=900, cmap='plasma')
Pxx, freqs, bins, im = ax2.specgram(x2, NFFT=NFFT, Fs=Fs, noverlap=900, cmap='plasma')
plt.show()








