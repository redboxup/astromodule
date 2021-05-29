#! /usr/bin/env python3.7

#program to plot the colour magnitude diagram 
#program to complete isochrone fitting on open cluster m67

from astromodule import *

#read the csv file into a dataframe
df = pd.read_csv('./output/cmd_data_m67.csv',
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp'])

def plot_cmd(df):
	y = df.phot_g_mean_mag.tolist()
	x = df.bp_rp.tolist()
	plt.plot(x,y, 'ro',markersize = 0.5,alpha =0.5)
	plt.ylabel ('$Absolute Magnitude(g)$')
	plt.xlabel ('$Color (bp-rp)$')
#	plt.xlim([xl,xm])
#	plt.ylim([yl,ym])	
	plt.gca().invert_yaxis()
	plt.title('Colour Magnitude Diagram')

#using plot_cmd function to plot the color magnitude diagram 
df0 = pd.read_csv("isochrone.csv")
color = df0.Gaia_BP - df0.Gaia_RP
plot_cmd(df)
#plt.plot( color.tolist(),df0.Gaia_BP.tolist())



#write a function to perform isochrone fitting on the above cmd
distance = 908
distance_modulus = 5*np.log10(distance/10)
#print ('dm = {}'.format(distance_modulus))
model_bp = df0.Gaia_BP.tolist()

model_bp = model_bp + distance_modulus
#plt.plot( color.tolist(),model_bp,'g')



#improving for extinction and reddening 
#ebv is reddening 
ebv = 0.041
#av = extinction 
av =3.2*ebv

model_bp = model_bp + av
#color = color.tolist()
color = color +ebv
plt.plot( color,model_bp)


plt.show()





