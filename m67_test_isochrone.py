#! /usr/bin/env python3.7

#program to test out isochrone fitting of a open cluster

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#importing darthmouth_isochrone group
from isochrones.dartmouth import Dartmouth_Isochrone

#now to read the data using pandas dataframe
df = pd.read_csv('./output/cmd_data_m67.csv', 
	usecols =['ra','dec','pmra','pmdec',
			'parallax','phot_g_mean_mag','bp_rp'])

#plot of the coluor magnitude diagram
fig , axis = plt.subplots(figsize=(8,6))
axis.plot(df.bp_rp,df.phot_g_mean_mag,'r.',alpha = 0.3)
axis.invert_yaxis()
axis.set_xlabel('bp_rp')
axis.set_ylabel('photo G mean magnitude')
plt.show()



