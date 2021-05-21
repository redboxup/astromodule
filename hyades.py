#! /usr/bin/env python3.7

#program to perform the astrometric analysis on the hyades star cluster

#read the csv file into a pandas dataframe
from astromodule import *
df = pd.read_csv("./database/hyades_3deg_photometric.csv",
	usecols=['ra','dec','pmra','pmdec',
	'parallax','phot_g_mean_mag','bp_rp','dr2_radial_velocity'])

#creating a dataframe to plot the propermotion of the stars in field of view
df1 = df.drop(columns=['ra','dec','parallax','phot_g_mean_mag','bp_rp',
			'dr2_radial_velocity'])
df1 = df1.dropna()

#To run gaussian mixture model on the data
def gmm_(df,n_components):
	X = df.to_numpy()
	gmm = mixture.GaussianMixture(n_components, covariance_type='full').fit(X)
	plot_results(X, gmm.predict(X), gmm.means_, gmm.covariances_,
		     'Gaussian Mixture')
	labels = gmm.predict(X)
	
	plt.show()
#	final_index = user_select(df,labels)
#	return final_index 	

#df2 = df1.loc[(df['pmra']>60)&(df['pmdec']<0)]

gmm_(df1,2)








