#! /usr/bin/env python3.7

#importing important query modules

from astroquery.gaia import Gaia
import seaborn as sns
from astropy.coordinates import SkyCoord
from astromodule import *

#now taking the galactic coordinates and converting it into equatorial 
l=[]
b = []
l0 = sys.argv[1]
b0 = sys.argv[2]
l.append(l0)
b.append(b0)

#now taking a box of values lying inside 0.1x0.1 deg square
coord = SkyCoord(l=l,b=b, unit=(u.degree,u.degree), frame='galactic')
coord = coord.icrs

#writing query in adql 
query = """SELECT ra,dec,pmra,pmdec,l,b,phot_g_mean_mag,parallax,dr2_radial_velocity,bp_rp FROM gaiaedr3.gaia_source 
           WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), 
           CIRCLE('ICRS',{ra},{dec},0.1))=1 ORDER BY random_index""".format(ra=str(coord.ra.deg[0]),dec=str(coord.dec.deg[0]))

#now launching the job
job = Gaia.launch_job_async(query)#,dump_to_file = True, output_format = 'csv')
r = job.get_results()
a = job.__dict__

#creating dataframe from a table
ralist = r['ra'].tolist()
declist = r['dec'].tolist()
pmralist = r['pmra'].tolist()
pmdeclist = r['pmdec'].tolist()
llist = r['l'].tolist()
blist = r['b'].tolist()
parallax_ =r['parallax'].tolist()
rv = r['dr2_radial_velocity'].tolist()
plist = r['phot_g_mean_mag'].tolist()
bprp  = r['bp_rp'].tolist()

df = pd.DataFrame(list(zip(pmralist,pmdeclist,llist,blist,ralist,declist,plist,parallax_,bprp,rv)),
                 columns = ['pmra','pmdec','l','b','ra','dec','G','parallax','bp_rp','dr2_radial_velocity'])
df1 = df

#Ending the job after creating the dataframe
job_id = a['jobid']
Gaia.remove_jobs(job_id)

#Because error for parallax is very high removing stars g>17
df1 = df1.drop(df1[df1.G >17].index)
no_stars = len(df1.G)

#Now df2 is the default starbase
df2 = df1
df1 = df1.drop(columns = ['dr2_radial_velocity'])
df1 = df1.dropna()

#printing total number of stars
print("Total number of stars are: {nox}".format(nox = no_stars))

plot1 = plt.figure(1)

pmra_ = df1.pmra.tolist()
pmdec_ = df1.pmdec.tolist()


plt.plot(pmra_,pmdec_,'ro',markersize = 1,alpha = 0.5)
plt.show(block = False)


# now in case we find the use of the 
ans = input("Do you need to implement density map? (y,n): ")
if ans == "n":
	print("density map will not be implemented")
else:
	#taking new limits
	xu = input("print upper limit for x: ")
	xl = input("print lower limit for x: ")
	yu = input("print upper limit for y: ")
	yl = input("print lower limit for y: ")
	xu = float(xu)
	xl = float(xl)
	yl = float(yl)
	yu = float(yu)

	#new dataframe becomes
	df1 = df1.drop(df1[df1.pmra < float(xl)].index)
	df1 = df1.drop(df1[df1.pmra > float(xu)].index)
	df1 = df1.drop(df1[df1.pmdec < float(yl)].index)
	df1 = df1.drop(df1[df1.pmdec > float(yu)].index)
	
	x = df1.pmra.to_numpy()
	y = df1.pmdec.to_numpy()
	plot2 = plt.figure(2)
	#plt.plot(pmra_,pmdec_,'ro',markersize = 1,alpha = 0.5)
	
	n_bins = 25
	h = plt.hist2d(x,y, bins = n_bins, cmap = 'Blues')
	cb = plt.colorbar()
	cb.set_label('count in bin')
	plt.show(block = False)
	n = input("how many points do you want to select: ")	
	pts = plt.ginput(n=int(n),show_clicks = True)
	
#	n_bins = float(n_bins)
	hx = (xu-xl)/n_bins
	hy = (yu-yl)/n_bins
	pts_listx = []
	pts_listy = []	
	for i in range(0,n_bins):
		pts_listx.append(xl+i*hx)
		pts_listy.append(yl+i*hy)
	
	#now we want to make find the nearest 
	def find_nearest(array,value):
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		return array[idx]
	val_box = []
	j = 0
	
	for i in pts:
		xbox = find_nearest(pts_listx,pts[j][0])
		ybox = find_nearest(pts_listy,pts[j][1])
		j = j+1
		val_box.append([xbox,ybox])
	
	#print(pts)
	#print(val_box)		

	#Now we can simpy use these to plot a color magnitude diagram
	#create a dataframe 
	j = 0
	for i in pts:
		df_temp = df1.drop(df1[df1.pmra < val_box[j][0]].index)
		df_temp = df1.drop(df1[df1.pmra > val_box[j][0]+hx].index)
		df_temp = df1.drop(df1[df1.pmdec < val_box[j][1]].index)
		df_temp = df1.drop(df1[df1.pmdec > val_box[j][1]+hy].index)
		plot = plt.figure(j+3)
		x= df_temp.bp_rp.tolist()
		y = df_temp.G.tolist()
		plt.plot(x,y, 'ro',markersize = 3,alpha =0.5)
		plt.ylabel ('$Absolute Magnitude(g)$')
		plt.xlabel ('$Color (bp-rp)$')
		plt.gca().invert_yaxis()
		plt.title('Colour Magnitude Diagram')
		plt.show(block = False)
		j = j+1
			

plt.show()































