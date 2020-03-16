 
"""
======================================
Spectrogram, power spectral density
======================================

Demo spectrogram and power spectral density on a frequency chirp.
"""

import numpy as np
import math
from matplotlib import pyplot as plt
import csv
from scipy import signal
from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale

"""
# Read file from slide
############################################################
with open('TR_2_slide_16.csv', 'r') as csv_file_slide:
    csv_reader_slide = csv.reader(csv_file_slide, delimiter=',')
    
    list_hx_s=[]
    list_hy_s=[]
    list_hz_s=[]
    line_count = 0
    
    for row in csv_reader_slide:
            list_hx_s.append(row[17])
            list_hy_s.append(row[18])
            list_hz_s.append(row[19])    
    
    del list_hx_s[0]
    del list_hy_s[0]
    del list_hz_s[0]

array_hx_s=np.array(list_hx_s, dtype='float64')
array_hy_s=np.array(list_hy_s, dtype='float64')
array_hz_s=np.array(list_hz_s, dtype='float64')

array_BT_s=np.sqrt(np.power(array_hx_s,2)+np.power(array_hy_s,2)+np.power(array_hz_s,2))
############################################################

# Generate a chirp signal
############################################################
sig_BT_s = array_BT_s
sig_BT2_s = np.sqrt(np.power((array_BT_s -array_BT_s.mean()),2)) #tHe variations


#plt.figure(figsize=(8, 5))
#plt.plot(sig_hx)

############################################################
# Compute and plot the spectrogram
############################################################
#
# The spectrum of the signal on consecutive time windows



freqs, times, spectrogram = signal.spectrogram(sig_BT_s)
plt.figure(figsize=(5, 4))
plt.imshow(spectrogram, aspect='auto', cmap='hot_r', origin='lower')
plt.title('Spectrogram')
plt.ylabel('Wavenumber band')
plt.xlabel('Space window')
plt.tight_layout()


############################################################
# Compute and plot the power spectral density (PSD)
############################################################
#
# The power of the signal per frequency band

freqs, psd = signal.welch(sig_BT2_s)

plt.figure(figsize=(5, 4))
plt.semilogx(freqs, psd)
plt.title('PSD: power spectral density')
plt.xlabel('Wavenumber')
plt.ylabel('Power')
plt.tight_layout()
plt.show()
############################################################



#Regression
############################################################
log_freqs=np.log(freqs)

#To remove nan and inf values
ii=np.isfinite(log_freqs)
log2=log_freqs[ii]
psd2=psd[1:] # Remove the first value 

slope, intercept, r_value, p_value, std_err = linregress(log2, psd2)
print(slope)
psd3=slope*log2 + intercept

"""
"""
fig1=plt.figure(1)
plt.semilogx(freqs, psd)
plt.semilogx(np.exp(log2), psd3, 'r')
plt.title('PSD: power spectral density')
#plt.text(0.05, 4, r'$\alpha_k='%slope)
plt.xlabel('Wavenumber ($kd_{i}$)')
plt.ylabel('Power')
plt.tight_layout()
plt.show(fig1)
"""
"""

fig1=plt.figure(1)
plt.semilogx(freqs, psd)
plt.semilogx(np.exp(log2), psd3, 'r')
plt.yscale('log')
plt.title('PSD: power spectral density')
#plt.text(60, .025, r'$\alpha_k='%slope)
plt.xlabel('Wavenumber ($kd_{i}$)')
plt.ylabel('Power')
plt.tight_layout()
plt.show(fig1)
############################################################
"""




# Read file from a box
############################################################

from matplotlib.backends.backend_pdf import PdfPages

#with open('TR_2_box_02.csv', 'r') as csv_file_box:
#with open('whis_box_80.csv', 'r') as csv_file_box:   
#with open('TR_3_lap_no_rand_phases.csv', 'r') as csv_file_box:       
#with open('TR_3_conty.csv', 'r') as csv_file_box:    
#with open('TR_3_conty_2.csv', 'r') as csv_file_box:
with open('TR_3_conty_2_100.csv', 'r') as csv_file_box:    
    csv_reader_box = csv.reader(csv_file_box, delimiter=',')
    
    list_hx_b=[]
    list_hy_b=[]
    list_hz_b=[]
    line_count = 0
    
    for row in csv_reader_box:
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
############################################################

# Generate a chirp signal
############################################################
sig_BT_b = array_BT_b 
sig_BT2_b = np.sqrt(np.power((array_BT_b -array_BT_b.mean()),2)) #The variations

wave_n_b, lenght_b, spectrogram_b = signal.spectrogram(sig_BT_b)


plt.figure(figsize=(5, 4))
plt.imshow(spectrogram_b, aspect='auto', cmap='hot_r', origin='lower')
plt.title('Magnetic field')
plt.ylabel('Wavenumber band')
plt.xlabel('Space window')
plt.tight_layout()

wave_n_2b, lenght_2b, spectrogram_2b = signal.spectrogram(sig_BT2_b)
plt.figure(figsize=(5, 4))
plt.imshow(spectrogram_2b, aspect='auto', cmap='hot_r', origin='lower')
plt.title('Variations of Bfield')
plt.ylabel('Wavenumber band')
plt.xlabel('Space window')
plt.tight_layout()


wave_n_b, psd_b = signal.welch(sig_BT_b)
plt.figure(figsize=(5, 4))
plt.semilogx(wave_n_b, psd_b)
plt.title('PSD: power spectral density')
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


"""
pp = PdfPages('foo.pdf')
pp.savefig(plot1)
pp.savefig(plot2)
pp.savefig(plot3)
pp.close()
"""

"""
fig, (ax1, ax2) = plt.subplots(nrows=2)
ax1.plot(t, x)
Pxx, freqs, bins, im = ax2.specgram(x, NFFT=NFFT, Fs=Fs, noverlap=900)
# The `specgram` method returns 4 objects. They are:
# - Pxx: the periodogram
# - freqs: the frequency vector
# - bins: the centers of the time bins
# - im: the matplotlib.image.AxesImage instance representing the data in the plot
plt.show()
"""

############################################################


    #a.pop(1)
    
#with open('whis_74.csv', mode='r') as csv_file:
#    csv_reader = csv.DictReader(csv_file)
#    for raw in csv_reader:
#        print(raw)

# Seed the random number generator
#np.random.seed(0)
#time_step = .01
#time_vec = np.arange(0, 70, time_step)
############################################################