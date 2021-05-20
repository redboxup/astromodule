#! /usr/bin/env python3.7

#program to perform astrometric analysis on the messier 67 open cluster

from astromodule import *
df = pd.read_csv("./database/m67_final.csv",
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp','dr2_radial_velocity'])

#creating a dataframe to plot the propermotion of the stars in field of view
df0 = df.drop(columns=
	['ra','dec','parallax','bp_rp','phot_g_mean_mag','dr2_radial_velocity'])

df0 = df0.dropna()

#running gmm and plot the results 
index = gmm_(df0,3)

#creating a dataframe to plot the vector plot diagrams 
df1 = df.drop(columns=
	['parallax','bp_rp','phot_g_mean_mag','dr2_radial_velocity'])
df1 = df1.dropna()

#filtering out the stars that belong to the cluster using indices
df2 = filter_cluster(df1,index)

#vector point diagram of the stars selected 
ax = plt.subplot()
ax.quiver(df2.ra, df2.dec,df2.pmra,df2.pmdec,scale=300)
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
print('data has been outputed to cmd_data_m67.csv')
df4.to_csv('./output/cmd_data_m67.csv') 

#now we output the data of galactic velocity co-ordinates 
df5 = df.dropna()
df5_index = df5.index.tolist()
mod_index = list(set(df5_index).intersection(index))
df5 = filter_cluster(df5,mod_index)
galactic_velocities(df5)







