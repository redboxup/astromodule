#! /usr/bin/env python3.7

#writing a python program to download the required databases and read from csv file

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.gaia import Gaia

#function to get square regions/areas

def square_reg(rax,decx,h,w):
	coord = SkyCoord(ra = rax, dec= decx, unit=(u.degree,u.degree),frame='icrs')
	width = u.Quantity(w,u.deg)
	height = u.Quantity(h,u.deg)
	Gaia.ROW_LIMIT = -1
	r = Gaia.query_object_async(coordinate = coord,width = width,height =height)
	return r

#calling the function to show the values

reg1 = square_reg(280,-60,0.1,0.1)
reg1.pprint()
