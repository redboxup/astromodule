#! /bin/bash
#opening the proper directory
cd ~/Desktop/git/astromodule/fortran

mkdir iso_temp

#run the fortran program to isolate the isochrone
./iso_interp_feh $1 $2 $3 $4 

#run the split fortran to split all the ages
cp isolf_split ./iso_temp/
cd ./iso_temp
./isolf_split main






