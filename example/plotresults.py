# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%% Dependencies

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%% Import results

filenames ={
    "u": "../output/u.txt",
    "xgrid": "../output/xgrid.txt"
    }

data = {}

for var, filename in filenames.items():
    data[var] = pd.read_csv(filename, delim_whitespace=True)
