#module for all the important functions that are being used in the project

#importing necessary modules
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
import subprocess as sb

from sklearn import mixture
from pandas import DataFrame
from scipy import linalg

#importing astropy modules
import astropy.coordinates as coord
import astropy.units as u

#declaring global variables
label_index = []
label_color = []
final_index = []



#function to plot the results of gmm and dpgmm
#function has been recursively called in the gmm_ and dpgmm_ function 
#ideally should not be used in the main file
def plot_results(X, Y_, means, covariances, title):
    color_iter = itertools.cycle(['b', 'g', 'r', 'c','y'])
    splot = plt.subplot(1, 1, 1)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        v, w = linalg.eigh(covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
        label_index.append(i)
        label_color.append(color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle)
        ell.set_clip_box(splot.bbox)
        ell.set_fill(0)
        splot.add_artist(ell)

    plt.title(title)

#Function to select the stars which belong to the open cluster
def user_select(df,labels):
	print(label_index,label_color)
	i=input("which label index is the cluster: ")
	i = int(i)
	for j in range(0,len(df)):
		if labels[j]==i:
			final_index.append(df.index[j])
	return final_index


#Fit a Gaussian mixture with EM using number of components
def gmm_(df,n_components):
	X = df.to_numpy()
	gmm = mixture.GaussianMixture(n_components, covariance_type='full').fit(X)
	plot_results(X, gmm.predict(X), gmm.means_, gmm.covariances_,
		     'Gaussian Mixture')
	labels = gmm.predict(X)
	
	plt.show()
	final_index = user_select(df,labels)
	return final_index 	




# Fit a Dirichlet process Gaussian mixture using five components
def dpgmm_(df,n_components):
	X = df.to_numpy()
	dpgmm = mixture.BayesianGaussianMixture(n_components,
		                                covariance_type='full').fit(X)
	plot_results(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_,
		     'Bayesian Gaussian Mixture with a Dirichlet process prior')
	labels = dpgmm.predict(X)
	
	plt.show()
	final_index = user_select(df,labels)
	return final_index 	



#function to plot vector plot diagram
def plot_vpd(df):
	ax = plt.subplot()
	ax.quiver(df.ra, df.dec,df.pmra,df.pmdec,300,color='k')
	plt.show()




#function to find the galactic velocity vectors
def galactic_velocities(df):
	df = df.dropna()
	l = len(df.ra)
	c = [0]*l
	gc = [0]*l
	print('galactic velocity coordinates have been written')
	sys.stdout = open("./output/galactic_vel_coordinates.dat","w")
	for i in range(0, l):
		j = df.index[i]
		c[i] = coord.SkyCoord(ra=df.ra[j]*u.degree, dec=df.dec[j]*u.degree,
		            distance=(df.parallax[j]*u.mas).to(u.pc, u.parallax()),
		            pm_ra_cosdec=df.pmra[j]*u.mas/u.yr,
		            pm_dec=df.dec[j]*u.mas/u.yr,
		            radial_velocity=df.dr2_radial_velocity[j]*u.km/u.s,
		            frame='icrs')

		gc[i] = c[i].transform_to(coord.Galactocentric)
		
		print(gc[i].v_x,gc[i].v_y,gc[i].v_z)
	sys.stdout.close()
	

#Now we will drop the rows that do not have the index in the list index
def filter_cluster(df,index):
	df_comp = df.drop(index)
	temp_arr= df_comp.index.tolist()
	df2 = df.drop(temp_arr)
	return df2





#using plot_cmd function to plot the color magnitude diagram 
def plot_cmd(df,xl,yl,xm,ym):
	y = df.phot_g_mean_mag.tolist()
	x = df.bp_rp.tolist()
	plt.plot(x,y, 'ro',markersize = 0.5,alpha =0.5)
	plt.ylabel ('$Absolute Magnitude(g)$')
	plt.xlabel ('$Color (bp-rp)$')
	plt.xlim([xl,xm])
	plt.ylim([yl,ym])	
	plt.gca().invert_yaxis()
	plt.title('Colour Magnitude Diagram')
	plt.show()



#obtain the relevant isochrone in a dataframe

def isochrone_(y,afe,feh,age):
	#running the bash script to run the fortran programs
	rc = sb.check_call(['./iso.sh', '16',str(y),str(afe),str(feh)])

	#read text file
	arr = os.listdir('./fortran/iso_temp')
	arr.sort()
	#print('The range of age is as follows')
	#print(arr[0],arr[-3])
	#print('\n')

	#now we create a file name
	def filename_gen(age):
		age = float(age)
		age = age*1000
		age = int(age)
		if age > 10000:
			age = str(age)
		else:
			age = str(age)
			age = '0' + age
		filename = 'a'+age+'main'
		return filename 

	file = filename_gen(age)
	if file in  arr:
		print('age is available')
	else:
		print('age unavailable')
		exit()

	#now we load the isochrone using numpy 
	fname = './fortran/iso_temp/'+file
	dtype_ = np.dtype([('EEP','i4'),('M/Mo', 'f4'),('logTeff','f4'),
		('LogG','f4'),('LogL/Lo','f4'),('Gaia_G','f4'),
		('Gaia_BP','f4'),('Gaia_RP','f4')])
	arr = np.loadtxt(fname ,dtype =dtype_,skiprows=9)

	#converting to dataframe 
	df = pd.DataFrame(arr)
	#print(df)

	#Extracting the data in headers
	with open(fname) as myfile:
	    head = [next(myfile) for x in range(9)]
	a=head[3]
	a = a.split()
	data = [[float(a[1]),float(a[2]),float(a[3]),
		float(a[4]),float(a[5]),float(a[6])]]
	cols = ['M_L','Y','Z','Zeff','FeH','aFe']
	df_h = pd.DataFrame(data,columns = cols)

	#removing the temporary file created during the process
	sb.check_call(["rm","-rf","./fortran/iso_temp"])
	return df,df_h





