#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:40:26 2019

@author: jaa
"""


import numpy as np
import math
from matplotlib import pyplot as plt
import csv
from scipy import signal
from scipy.stats import linregress
from matplotlib.ticker import NullFormatter  # useful for `logit` scale


number_cores=[16,32,64,128,256,512]
walltime_sec=[247119, 119674, 56770, 30580, 14231.9, 7867.95]

walltime_hours=[x/3600 for x in walltime_sec]
walltime_hours2=[x/(1/512) for x in walltime_hours]

idealtime=[1/x for x in number_cores]
idealtime2=[x/(1/1024) for x in idealtime] #why 1024 fix?


f=plt.figure(figsize=(5, 4))

f=plt.figure()
plt.plot(number_cores, idealtime, '-+')
plt.plot(number_cores, walltime_hours, '-*')
#plt.title('weak scalling')
plt.ylabel('Wall time (Hours)')
plt.xlabel('Number of cores')
plt.tight_layout()
plt.show()


f2, ax = plt.subplots()
plt.loglog(number_cores, idealtime2, '-')
plt.loglog(number_cores, walltime_hours, '-*')
#plt.title('weak scalling')
plt.ylabel('Wall time (Hours)')
plt.xlabel('Number of cores')
#ax.set_yscale('log', basex=2)
ax.set_xscale('log', basex=2)
plt.gca().legend(('Ideal Scaling','PSC Scaling'))
plt.tight_layout()
plt.show()



f.savefig("walltimevscores.png", bbox_inches='tight')
f2.savefig("walltimevscores_loglog_2.png", bbox_inches='tight')