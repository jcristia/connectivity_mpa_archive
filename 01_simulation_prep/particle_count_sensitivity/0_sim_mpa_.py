# Opendrift simulations to test the sensitivity of the amount of particles
# released. I am following the approach from Simons et al. 2013.
# More detailed notes in the scripts in the analysis folder.


sims = 20
mpa_uid = '207'
dates = '20110801_20111014'


from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_shape
from datetime import datetime
from datetime import timedelta
from osgeo import ogr
import numpy as np
import os

######## Set up directory structure

base = os.path.dirname(__file__)
outputs = os.path.join(base, 'outputs')
nc_out = os.path.join(outputs, 'nc')
np_out = os.path.join(outputs, 'npy')
logs = os.path.join(outputs, 'logs')
imgs = os.path.join(outputs, 'img')


######## Readers

nc_nep = 'NEP36_1h_{}.nc'.format(dates)
nc_ssc = 'SalishSea_1h_{}_opendrift.nc'.format(dates)

base_mod = 'C:/Users/jcristia/Documents/GIS/MSc_Projects/Hakai/spatial/models'
base_lnd = 'inputs'
shp_lnd = 'landmask_FINAL_wgs84.shp'

file_nep = os.path.join(base_mod, 'nep_nemo/processed', nc_nep)
file_ssc = os.path.join(base_mod, 'salishsea/salishseacast/forcing', nc_ssc)
file_lnd = os.path.join(base, base_lnd, shp_lnd)
reader_nep = reader_netCDF_CF_unstructured.Reader(
    file_nep, latstep=0.01, lonstep=0.01, buffer=0.1, name='NEP')
reader_ssc = reader_netCDF_CF_unstructured.Reader(
    file_ssc, latstep=0.004, lonstep=0.004, buffer=0.1, name='SSC')
reader_lnd = reader_shape.Reader.from_shpfiles(file_lnd)
# buffers are set low and have been tested to not get data block errors
# for boundaries and around certain MPAs you will still get warnings about
# extrapolation. This is ok. See main note in Evernote for details.


######## Get MPA and particle number

mpa_shp = os.path.join(base, 'inputs/mpa_.shp')
shp = ogr.Open(mpa_shp)
lyr = shp.GetLayer(0)
for feature in lyr:
    uID = feature.GetField('uID_202011')
    if uID == int(mpa_uid):
        featurenum = feature.GetFID() + 1 # opendrift subtracts 1 for some reason
        uID_202011 = feature.GetField('uID_202011')
        particles = feature.GetField('part_num_T')
print('Feature number: ' + str(featurenum))
print('Particles: ' + str(particles))


######## Iterate through simulations

for sim in range(sims):

    ######## Model

    o = OceanDrift(
        loglevel=20, 
        logfile= os.path.join(logs, 'log_{}_{}.log'.format(mpa_uid, str(sim))),
        seed=None  # Even though I am doing sensitivity tests, this still has to
        # be set at None. I controls starting positions (which would be fine if 
        # they were constant), but it also controls diffusion, which I need to
        # be random or else everything will be the same and defeat the point.
        )
    o.add_reader([reader_lnd, reader_ssc, reader_nep])


    ######### Seed particles

    time_step = timedelta(hours=4)
    num_steps = 84
    for i in range(num_steps):
        o.seed_from_shapefile(
            mpa_shp, # this didn't work if I did shp here instead of mpa_shp
            featurenum = featurenum,
            origin_marker=uID_202011,
            number=particles, 
            time=reader_ssc.start_time + i*time_step
            )


    ######### starting coordinates, for use in biology script

    if sim<2: # just do this for the first 2 to see how much they differ
        npy_lon = os.path.join(np_out, 'lon_{}_{}.npy'.format(mpa_uid, str(sim)))
        npy_lat = os.path.join(np_out, 'lat_{}_{}.npy'.format(mpa_uid, str(sim)))    
        np.save(npy_lon, o.elements_scheduled.lon)
        np.save(npy_lat, o.elements_scheduled.lat)


    ######### Configure and Run

    #o.list_configspec()
    o.set_config('general:use_auto_landmask', False)  # use custom landmask
    o.set_config('drift:current_uncertainty', 0) # using basemodel hardcoded values
    o.set_config('general:coastline_action', 'stranding')
    o.set_config('drift:scheme', 'runge-kutta')

    output_nc = os.path.join(nc_out, 'output_{}_{}.nc'.format(mpa_uid, str(sim)))
    print("Simulation for sim {}".format(str(sim)))
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

    o.plot(filename=os.path.join(imgs, 'output_{}_{}.jpg'.format(mpa_uid, str(sim))))

    #####################