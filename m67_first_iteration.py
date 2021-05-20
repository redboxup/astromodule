#! /usr/bin/env python3.7

#program to perform complete astrometric analysis on the m67 open cluster

#importing necessary libraries
from astromodule import *

#reading the csv file
d = pd.read_csv('./database/m67_gaia_edr3_6179rows.csv',usecols=['pmra','pmdec'])
d = d.dropna()
d = d.reset_index(drop=True)

#calling the gaussian mixture modelling function
plot1 = plt.figure(1)
gmm_(d,n_components=3)

#plotting the vector plot diagram of the selected cluster of stars
#temporary we are using pm_check file
df = pd.read_csv('./database/m67_pm_check.csv',usecols =['ra','dec','pmra','pmdec','dr2_radial_velocity','parallax'])
plot2 = plt.figure(2)
ax = plt.subplot()
ax.quiver(df.ra, df.dec,df.pmra,df.pmdec,scale=300)
plt.show()
	



#calculate the galactic velocities
galactic_velocities(df)

