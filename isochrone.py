#! /usr/bin/env python3.7

#program to plot the isochrones 
from astromodule import *

y = 1
afe = 2
feh = 0.041
age = 1.5
df,df_h = isochrone_(y,afe,feh,age)
