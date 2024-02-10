# compile the nep36 netcdf files based on a date range


import netCDF4 as nc
import numpy as np
import xarray as xr
import os
import datetime
import cftime
from dateutil.parser import parse


#######
# Select files by date
#######

year = '2014'
# from and to (enter month and days as 2 digits)
fr = '0101'
to = '0316'

# get the u and v files that contain these dates
paths = {
    'base': r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\raw',
    'mesh': r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\mesh_mask.nc',
    'out': r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\processed'
}
fri = int(year + fr)
toi = int(year + to)
ufiles = []
vfiles = []
for vel in ['U', 'V']:
    for file in os.listdir(os.path.join(paths['base'],year)):
        #print(file)
        uv = file.split('_')[5]
        if uv == vel:
            d_range = file.split('_')[7]
            f = int(d_range.split('-')[0])
            t = int((d_range.split('-')[1]).split('.')[0])
            if fri >= f and fri <= t:
                if uv == 'U':
                    ufiles.append(file)
                else:
                    vfiles.append(file)
            elif toi >= t:
                if uv == 'U':
                    ufiles.append(file)
                else:
                    vfiles.append(file)
            elif toi >= f: # this should catch the last one
                if uv == 'U':
                    ufiles.append(file)
                else:
                    vfiles.append(file)


#######
# Format
#######

# Load forcing data
u = []
for file in ufiles:
    with xr.open_dataset(os.path.join(paths['base'], year, file)) as data:
        u.append(data.uos.values)
u = np.array(u)
u = u.reshape(-1, *u.shape[2:]) # this flattens the first two dimensions so that it is no longer divided out by file, just by time
v = []
for file in vfiles:
    with xr.open_dataset(os.path.join(paths['base'], year, file)) as data:
        v.append(data.vos.values)
v = np.array(v)
v = v.reshape(-1, *v.shape[2:])


# Daterange for simulation
dates = [year + fr + ' ' + '00:30', year + to + ' ' + '23:30']
daterange = [parse(d) for d in dates]
times = []
for file in ufiles: # only need to get time once
    with xr.open_dataset(os.path.join(paths['base'], year, file)) as data:
        times.append(data.time_counter.values)
times = np.array(times)
times = times.reshape(-1) # flatten it so that it is not divided into different arrays by file
timesds = xr.Dataset(coords = {'time': times})  # put this back as an xarray dataset so that I can use 'sel'
time = timesds.time.sel(time=slice(*daterange))
# 'time' is the actual dates I will use as the dimension
# but to select the uv I want by time I also need a boolean of indexes to select:
fnp = np.datetime64(daterange[0])
tnp = np.datetime64(daterange[1])
tx = np.where((times >= fnp) & (times <= tnp), True, False)


# lat/lon
grid = xr.open_dataset(os.path.join(paths['base'], year, ufiles[0]))


# land mask
mask = xr.open_dataset(os.path.join(paths['mesh']))
tmask = mask.tmask[0, 0, :, :].values.reshape(-1).astype(bool) # just need the first time step and surface layer


# out path
fn = 'NEP36_1h_' + '_'.join(d.strftime('%Y%m%d') for d in daterange) + '.nc'
forcing_path = os.path.join(paths['out'], fn)

# reshape - select time slice, remove landpoints
ush = u.reshape(np.shape(u)[0], -1)
ush = ush[tx]
ush = ush[:,tmask]
vsh = v.reshape(np.shape(v)[0], -1)
vsh = vsh[tx]
vsh = vsh[:,tmask]

# Reshape, remove landpoints, and save to local netCDF path
ds = xr.Dataset(
    {
        'longitude': ('flat', grid.nav_lon.values[:].reshape(-1)[tmask]), # the flat references that it is the total amount of lat-lon pairs that the u and v are represented along, as opposed to separately. It is the dimension amount.
        'latitude': ('flat', grid.nav_lat[:].values.reshape(-1)[tmask]),
        'u': (['time', 'flat'], ush, {'standard_name': 'x_sea_water_velocity'}),
        'v': (['time', 'flat'], vsh, {'standard_name': 'y_sea_water_velocity'}),
    },
    coords={'time': time}
    # xarray reformats time so that the origin is just the start of the dataset and the units are the smallest difference in time step (which is hours in this case). Very handy.
).to_netcdf(forcing_path)


# manual check (Arc can't read it directly when lat/lon aren't dimensions)
# import numpy as np
# import xarray as xr
# filename = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\processed\NEP36_1h_20110101_20110316.nc"
# chf = xr.open_dataset(filename)
# lat = chf.latitude[:].values
# lon = chf.longitude[:].values
# uo = chf.u[0].values
# vo = chf.v[0].values
# headers = ['X', 'Y', 'u', 'v']
# output = open(r'C:\Users\jcristia\Desktop\grid_nepnemo_JC.csv', 'w')
# for header in headers:
#     output.write(header + ",")
# output.write("\n")
# for x,y,u,v in zip(lon,lat,uo,vo):
#     if x != 0 and y != 0:
#         # can only write strings
#         output.write('%s' % x)
#         output.write(',')
#         output.write('%s' % y)
#         output.write(',')
#         output.write('%s' % u)
#         output.write(',')
#         output.write('%s' % v)
#         output.write(',')
#         output.write("\n")
# output.close()