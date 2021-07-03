#! /usr/bin/env python3.7

#importing important query modules

from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import os

#now taking the galactic coordinates and converting it into equatorial 
l=[]
l0 = input("galactic longitude: ")
l.append(l0)
b = [0.5]

#now taking a box of values lying inside 0.1x0.1 deg square
coord = SkyCoord(l=l,b=b, unit=(u.degree,u.degree), frame='galactic')
coord = coord.icrs

#writing query in adql
query = """SELECT ra,dec,pmra,pmdec FROM gaiaedr3.gaia_source 
           WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), 
           BOX('ICRS',{ra},{dec},0.1,0.1))=1 ORDER BY random_index""".format(ra=str(coord.ra.deg[0]),dec=str(coord.dec.deg[0]))

#going to the proper directory where we need to save the csv file
os.chdir('./skydb')

#now launching the job
job = Gaia.launch_job_async(query,dump_to_file = True, output_format = 'csv')
r = job.get_results()
a = job.__dict__
oldfilename = a['outputFile']
newfilename = 'long'+str(l[0])+'_lat'+str(b[0])+'.csv'

#renaming the file 
os.rename(oldfilename,newfilename)

#going back to other directory to save the plots
os.chdir('..')
os.chdir('./skyplot')

#plotting the result
ralist = r['ra'].tolist()
declist = r['dec'].tolist()
pmralist = r['pmra'].tolist()
pmdeclist = r['pmdec'].tolist()

#creating a dataframe to store the values
df = pd.DataFrame(list(zip(ralist,declist,pmralist,pmdeclist)),
                 columns = ['ra','dec','pmra','pmdec'])
file_ra_dec = 'long'+str(l[0])+'_lat'+str(b[0])+'_ra_dec.png'
filepm_ra_dec = 'long'+str(l[0])+'_lat'+str(b[0])+'_pmra_pmdec.png'

seaborn.jointplot(x='ra',y='dec',data=df,s = 3,alpha = 0.5)
plt.savefig(file_ra_dec,dpi = 500)

seaborn.jointplot(x='pmra',y='pmdec',data = df,s = 3,alpha = 0.5)
plt.savefig(filepm_ra_dec,dpi = 500)
plt.show()


