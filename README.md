# Summary of work: Astrometric Analysis of Open Clusters
 
----------------------------------------------

## Abstract

The primary objective of the project was to study and identify the astrometric paprameters of M67 and Hyades open cluster, using data provided by ESO Gaia Survey in Early Data Relase 3 (EDR3). Firstly, most probable members of the respective open clusters were identified by proper motion analysis. Next, the Colour Magnitude Diagram (CMD) of the cluster is created. Then a theoretical isochrone curve is fitted to determine important astrometric parameters of the cluster.  


## Methodology

### Collection of Data
Gaia is a European space mission providing astrometry, photometry and spectroscopy for more than 1 billion stars in the Milky Way. We can get the information of all the stars lying within radius of some arc seconds with center at the specified right ascension (RA) and declination (DEC). Since M67 and Hyades are well known open clusters the RA and DEC of centers of open clusters were known. 

## Proper Motion Analysis

We first performed a proper motion analysis. We plot the pmRA along the horizontal axis and pmDEC along the vertical axis. We see a high density clump in the figure. This is indicative of presence of open cluster. We use Gaussian Mean Mixture (GMM) model to select all the stars which most probably belong to the open cluster.

<p align = "center">
<img src="https://github.com/redboxup/astromodule/blob/main/plots/hyades_gmm_plot.png" width="425"/> <img src="https://github.com/redboxup/astromodule/blob/main/plots/m67_gmm_plot.png" width="415"/> 
</p>

We can verify that we have members from the open cluster we can plot the stars proper motion as velocities from their positions with RA and DEC as horizontal and veritical axis repectively.

<p align = "center">
<img src="https://github.com/redboxup/astromodule/blob/main/plots/pm_plot_hyades.png" width="395"/> <img src="https://github.com/redboxup/astromodule/blob/main/plots/pm_m67.png" width="425"/> 
</p>






