#! /usr/bin/env python3.7

#writing a python program to download the required databases and read from csv file

from astroquery.gaia import Gaia

job = Gaia.launch_job_async("select top 100 designation,ra,dec "
                             "from gaiadr2.gaia_source order by source_id")
r = job.get_results()
print(r)
