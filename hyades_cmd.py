#! /usr/bin/env python3.7

#program to plot the colour magnitude diagram 
#program to complete isochrone fitting on open cluster m67

from astromodule import *

#read the csv file into a dataframe
df = pd.read_csv('./output/cmd_data_hyades.csv',
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp'])

def plot_cmd(df):
	y = df.phot_g_mean_mag.tolist()
	x = df.bp_rp.tolist()
	plt.plot(x,y, 'ro',markersize = 1,alpha =0.5)
	plt.ylabel ('$Absolute Magnitude(g)$')
	plt.xlabel ('$Color (bp-rp)$')
#	plt.xlim([xl,xm])
#	plt.ylim([yl,ym])	
	plt.gca().invert_yaxis()
	plt.title('Isochrone fitting')

#writing a function that will be repeated mutiple times


def isochrone_plot(y,afe,feh,ebv,age,distance):
	#using plot_cmd function to plot the color magnitude diagram 
	df0,df0_h = isochrone_(y,afe,feh,age)
	color = df0.Gaia_BP - df0.Gaia_RP
	plot_cmd(df)

	#write a function to perform isochrone fitting on the above cmd
	distance_modulus = 5*np.log10(distance/10)
	print ('dm = {}'.format(distance_modulus))
	print ('age = {}'.format(age))
	print(df0_h)
	model_bp = df0.Gaia_BP.tolist()
	model_g  = df0.Gaia_G.tolist()
	model_bp = model_bp + distance_modulus
	model_g = model_g + distance_modulus


	#taking into consideration extinction and reddening 
	#ebv is reddening 
#	ebv = 0.04
	#av is extinction 
	av =3.2*ebv

	model_g = model_g +av
	model_bp = model_bp + av
	color = color +ebv
	plt.plot( color,model_g)


y = 1
afe = 2
feh = sys.argv[1]	#feh = 0
ebv = sys.argv[2]	#ebv = 0.04	
age= sys.argv[3]	#age = 4
distance = sys.argv[4]	#distance = 851
isochrone_plot(y,afe,float(feh),float(ebv),float(age),float(distance))
#plt.pause(5)
plt.show()
#var = 1
#while var ==1:
	



