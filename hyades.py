#! /usr/bin/env python3.7

#program to perform the astrometric analysis on the hyades star cluster

#read the csv file into a pandas dataframe
from astromodule import *
from astropy.coordinates import Distance

df = pd.read_csv("./database/hyades_3deg_photometric.csv",
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp','dr2_radial_velocity'])
#creating a dataframe to plot the propermotion of the stars in field of view
df1 = df.drop(columns=['phot_g_mean_mag','bp_rp',
			'dr2_radial_velocity'])
df1 = df1.dropna()

#removing from the cluster field of view stars 

#first to find the distance of all the stars
df1_parallax0 = df1.parallax.to_list()
df1_parallax1 = list(map(abs,df1_parallax0))
d = Distance(parallax = df1_parallax1*u.mas)
d = d.to_value()

#adding the distance column to the dataframe
df1['distance']=d

#now removing all the rows with a distance >100 pc
df1 = df1.drop(df1[df1.distance >100].index)
df2 = df1.drop(columns = ['parallax','distance','ra','dec'])

#running the gmm algorithm
index = gmm_(df2,2)

#filtering and using the stars that were selected
df3 = filter_cluster(df1,index)


#plotting the vector plot diagram for the given star cluster
ax = plt.subplot()
ax.quiver(df3.ra, df3.dec,df3.pmra,df3.pmdec,scale = 1000,color='k')
plt.show()

#creating a dataframe for creating a color magnitude diagram
df3 =df.drop(columns=
	['dr2_radial_velocity']) 
df3 = df3.dropna()
df3_index = df3.index.tolist()

#removing indices with no corresponding radial velocites
mod_index = list(set(df3_index).intersection(index))

#filtering out the stars that belong to the cluster using indices 
df4 = filter_cluster(df3,mod_index)

#now we output the data to a csv file
print('data has been outputed to cmd_data_hyades.csv')
df4.to_csv('./output/cmd_data_hyades.csv') 





