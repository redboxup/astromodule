#! /usr/bin/env python3.7

#program to plot the colour magnitude diagram 
#program to complete isochrone fitting on open cluster m67

from astromodule import *

#read the csv file into a dataframe
df = pd.read_csv('./output/cmd_data_m67.csv',
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp'])


#using plot_cmd function to plot the color magnitude diagram 
plot_cmd(df,0,8,3,22)


#write a function to perform isochrone fitting on the above cmd


