#! /usr/bin/env python3.7

#program to test out how to spit labels to clusters of information

from astromodule import *

#creating a dataframe to plot the propermotion of the stars in field of view
df0 = pd.read_csv("./database/m67_final.csv",usecols =['pmra','pmdec'])
df0 = df0.dropna()

#running gmm and plot the results 
index = gmm_(df0,3)

#creating a dataframe to plot the vector plot diagrams 
df1 = pd.read_csv('./database/m67_final.csv',usecols =['ra','dec','pmra','pmdec'])
df1 = df1.dropna()

#filtering out the stars that belong to the cluster using indices
df2 = filter_cluster(df1,index)

#vector point diagram of the stars selected 
ax = plt.subplot()
ax.quiver(df2.ra, df2.dec,df2.pmra,df2.pmdec,scale=300)
plt.show()


#creating a dataframe for creating a color magnitude diagram
df3 = pd.read_csv("./database/m67_final.csv",usecols=['ra','dec','pmra','pmdec','parallax','phot_g_mean_mag','bp_rp'])
df3 = df3.dropna()
df3_index = df3.index.tolist()

#removing indices with no corresponding radial velocites
mod_index = list(set(df3_index).intersection(index))

#filtering out the stars that belong to the cluster using indices 
df4 = filter_cluster(df3,mod_index)

#now we output the data to a csv file
print('data has been outputed to cmd_data_m67.csv')
df4.to_csv('./output/cmd_data_m67.csv') 











