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
print(df1)


#plotting the pmra and pmdec without performing gmm
def plot_pm(df):
	y = df.pmdec.tolist()
	x = df.pmra.tolist()
	plt.plot(x,y, 'ro',markersize = 0.5,alpha =0.5)
	plt.show()
#plot_pm(df1)

df2 = df1.drop(columns = ['parallax','distance','ra','dec'])

#To run gaussian mixture model on the data
def gmm_(df,n_components):
	X = df.to_numpy()
	gmm = mixture.GaussianMixture(n_components, covariance_type='full').fit(X)
	plot_results(X, gmm.predict(X), gmm.means_, gmm.covariances_,
		     'Gaussian Mixture')
	labels = gmm.predict(X)
	
	plt.show()
	final_index = user_select(df,labels)
	return final_index 	


index = gmm_(df2,2)

df3 = filter_cluster(df1,index)

ax = plt.subplot()
ax.quiver(df3.ra, df3.dec,df3.pmra,df3.pmdec,scale = 1000,color='k')
plt.show()







