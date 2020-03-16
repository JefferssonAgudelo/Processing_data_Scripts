#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:31:06 2019

@author: jaa
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import viscid
import pandas as pd
viscid.calculator.evaluator.enabled = True
from viscid.plot import vpyplot as vlt



# This load the variables
import h5_filereader_jeff as reader
reader.ex.shape



run = "/run/media/jaa/C2BCB9BCBCB9AB75/DIRAC_RUNs/SCALING/weak_2_scaling_256/pfd.xdmf"
vf = viscid.load_file(run, force_reload=True)