#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 21:39:55 2019

@author: oem
"""

import h5py
import numpy as np

filename = '/home/oem/PIC_STUFF/Whistler_1/pfd.000200_p000000.h5'
f = h5py.File(filename, 'r')

# List all groups

# There are 25 groups in each file.h5
list (f.keys())

for name in f: print(name) # This print the names (as list)

print("Keys: %s" % f.keys())
Temperature = list(f.keys())[0] #This is the temperture data
B_field = list(f.keys())[17] #This is the magnetic field data

# Get the data
list_Tempe = list(f[Temperature])
list_B_field = list(f[B_field])




# How to write

data_matrix = np.random.uniform(-1, 1, size=(10, 3))

# Write data to HDF5
data_file = h5py.File('file.hdf5', 'w')
data_file.create_dataset('group_name', data=data_matrix)
data_file.close()
 
