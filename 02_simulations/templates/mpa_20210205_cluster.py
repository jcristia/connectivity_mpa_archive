


from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_shape
from datetime import datetime
from datetime import timedelta
from osgeo import ogr
import numpy as np
import os

######## Set up directory structure

# when running on Cedar in a singularity image the cwd gets changed back to
# /home/jcristia (or whatever I set the singularity exec --home tag to)
# I'm just going to build file paths from root
root = os.getcwd()
base = os.path.dirname(__file__) # on Cedar, this will give me the path relative to where I execute the job from or where I set the singularity exec --home tag to
outputs = os.path.join(root, base, 'outputs')
nc_out = os.path.join(outputs, 'nc')
np_out = os.path.join(outputs, 'npy')
logs = os.path.join(outputs, 'logs')
for d in [outputs, nc_out, np_out, logs]:
    if not os.path.exists(d):
        os.makedirs(d)


######## Model

o = OceanDrift(
    loglevel=20, 
    logfile= os.path.join(logs, 'log_{}.log'.format(shp_group)),
    seed=None
    )


######## Readers

nc_nep = 'NEP36_1h_{}.nc'.format(dates)
nc_ssc = 'SalishSea_1h_{}_opendrift.nc'.format(dates)

base_mod = 'inputs/hydro_models'
base_lnd = 'inputs/coastline'
shp_lnd = 'landmask_FINAL_wgs84.shp'

file_nep = os.path.join(root, base_mod, nc_nep)
file_ssc = os.path.join(root, base_mod, nc_ssc)
file_lnd = os.path.join(root, base_lnd, shp_lnd)
reader_nep = reader_netCDF_CF_unstructured.Reader(
    file_nep, latstep=0.01, lonstep=0.01, buffer=0.1, name='NEP')
reader_ssc = reader_netCDF_CF_unstructured.Reader(
    file_ssc, latstep=0.004, lonstep=0.004, buffer=0.1, name='SSC')
reader_lnd = reader_shape.Reader.from_shpfiles(file_lnd)
# buffers are set low and have been tested to not get data block errors
# for boundaries and around certain MPAs you will still get warnings about
# extrapolation. This is ok. See main note in Evernote for details.

o.add_reader([reader_lnd, reader_ssc, reader_nep])


######### Seed particles

# for loops must be in this order so that the particles are placed in order by
# when they are seeded, not in order by feature

time_step = timedelta(hours=4)
num_steps = 84
for i in range(num_steps):
    mpa_shp = os.path.join(root, 'inputs/mpas/mpa_.shp')
    shp = ogr.Open(mpa_shp)
    lyr = shp.GetLayer(0)
    for feature in lyr:
        group = feature.GetField('grouping')
        if group == shp_group:

            featurenum = feature.GetFID() + 1 # opendrift subtracts 1 for some reason
            uID = feature.GetField('uID_202011')
            particles = feature.GetField('part_num')
            part_fact = 1 # factor to reduce particle count for testing
            particles = int(particles * part_fact)

            o.seed_from_shapefile(
                mpa_shp, # this didn't work if I did shp here instead of mpa_shp
                featurenum = featurenum,
                origin_marker=uID,
                number=particles, 
                time=reader_ssc.start_time + i*time_step
                )

# starting coordinates, for use in biology script
npy_lon = os.path.join(np_out, 'lon_{}.npy'.format(shp_group))
npy_lat = os.path.join(np_out, 'lat_{}.npy'.format(shp_group))    
np.save(npy_lon, o.elements_scheduled.lon)
np.save(npy_lat, o.elements_scheduled.lat)


######### Configure and Run

#o.list_configspec()
o.set_config('general:use_auto_landmask', False)  # use custom landmask
o.set_config('drift:current_uncertainty', 0) # using basemodel hardcoded values
o.set_config('general:coastline_action', 'stranding')
o.set_config('drift:scheme', 'runge-kutta')

output_nc = os.path.join(nc_out, 'output_{}.nc'.format(shp_group))
print("Simulation for group {}".format(shp_group))
o.run(
    #steps=720,
    end_time=reader_nep.end_time, 
    time_step=120, 
    time_step_output=1800, 
    outfile= output_nc, 
    export_variables=[
        'age_seconds', 
        'land_binary_mask',
        'origin_marker']
    )
print(o)


######### Outputs

#o.plot(filename='output_.jpg')
#o.animation()

#####################