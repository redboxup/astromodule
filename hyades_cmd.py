#! /usr/bin/env python3.7

#program to plot the colour magnitude diagram 
#program to complete isochrone fitting on open cluster hyades

from astromodule import *

#read the csv file into a dataframe
df = pd.read_csv('./output/cmd_data_hyades.csv',
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp'])


#using plot_cmd function to plot the color magnitude diagram 
def plot_cmd(df):
	y = df.phot_g_mean_mag.tolist()
	x = df.bp_rp.tolist()
	plt.plot(x,y, 'ro',markersize = 1,alpha =0.5)
	plt.ylabel ('$Absolute Magnitude(g)$')
	plt.xlabel ('$Color (bp-rp)$')
	#plt.xlim([xl,xm])
	#plt.ylim([yl,ym])	
	plt.gca().invert_yaxis()
	plt.title('Colour Magnitude Diagram')
	plt.show()

plot_cmd(df)


#write a function to perform isochrone fitting on the above cmd

