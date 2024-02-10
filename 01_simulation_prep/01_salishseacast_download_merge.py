
# Purpose:
# download and merge salishseacast data
# there is a 2GB limit on downloads
# so need to download smaller chunks and merge

import xarray as xr
import os
from datetime import datetime, timedelta
from calendar import monthrange
import numpy as np
from tqdm import tqdm


##########
# helper functions
##########

def load_paths():
    paths = {
        'erddap': 'https://salishsea.eos.ubc.ca/erddap/griddap',
        'out': r'D:\Hakai\models\salishsea\salishseacast\ssc_',
    }
    return paths

def load_netCDF_keys(filesystem='errdap'):
    # NetCDF file keys master dict
    if filesystem == 'errdap':
        key_dict = {
            'temperature': 'g3DTracer',
            'salinity': 'g3DTracer',
            'nitrate': 'g3DBiology',
            'uVelocity': 'g3DuGrid',
            'vVelocity': 'g3DvGrid',
            'u_wind': 'aSurfaceAtmosphere',
            'v_wind': 'aSurfaceAtmosphere',
        }        
    return key_dict

def extract_variables(data_vars, ds, variables, key_dict, dates=[None], dim='time', indices={'gridX': slice(None), 'gridY': slice(None)}):

    # Define time index dict
    tindex = {dim: slice(*dates)}  # * unpacks list items, slice creates indices out of those dates.
    
    for var in variables:
        # Initialize data array
        if not var in data_vars:
            data_vars[var] = ds[key_dict[var]][var].isel(indices).sel(tindex).load()
        else: # Concatenate data arrays
            data_vars[var] = xr.concat([data_vars[var], ds[key_dict[var]][var].isel(indices).sel(tindex).load()], dim=dim,)        
    return data_vars


########
# Main function
########

def extract_hindcast(daterange, variables, res='1h', version='19-05', filesystem='errdap', indices={'x': slice(None), 'y': slice(None)}):
    
    # Prepare variable definitions
    years, months, days = [[getattr(date, key) for date in daterange] for key in ['year', 'month', 'day']]
    paths = load_paths()
    key_dict = load_netCDF_keys(filesystem=filesystem)
    ds, keys = {}, list(set([key_dict[var] for var in variables]))
    encoding = dict(zip(variables, np.repeat({'zlib': True}, len(variables))))
    prefix_out = os.path.join(paths['out'])
    
    # Initiate loading protocol based on filesystem
    if filesystem == 'errdap':
        #prefix_out = f'{prefix_out}{key_dict[variables[0]][1:]}_'
        for key in keys:
            ds[key] = xr.open_dataset(paths['erddap'] + f'/ubcSS{key}Fields{res}V{version}')
        attrs = ds[key_dict[variables[0]]].attrs
    else:
        raise ValueError(f'Unknown filesystem: {filesystem}')

    # Loop through years
    for year in range(years[0], years[1] + 1):
        
        # Initialize data_vars dict and parse months
        data_vars = {}
        monthday = [[1, 1], [12, 31]]
        monthspan = [1, 13]
        if year == years[0]: monthspan[0] = months[0]
        if year == years[1]: monthspan[1] = months[1] + 1
            
        # Extract data month by month
        for month in tqdm(range(*monthspan), desc=f'Loading {year}'):

            # Parse daterange
            day, monthdays = 1, monthrange(year, month)[1]
            if (year == years[0]) and (month == months[0]):
                day = days[0]
                monthdays = monthdays - day + 1
                monthday[0] = [month, day]
            if (year == years[1]) and (month == months[1]):
                monthdays = days[1]
                monthday[1] = [month, monthdays]
            startdate = datetime(year, month, day)

            # Load variables from ERDDAP using specified month range
            if filesystem == 'errdap':
                dates = [startdate, startdate + timedelta(monthdays)]
                data_vars = extract_variables(data_vars, ds, variables, key_dict, dates=dates, indices=indices)
            else:
                raise ValueError(f'Unknown filesystem: {filesystem}')

        # Save as netCDF file
        datestr = '_'.join(datetime(year, *md).strftime('%Y%m%d') for md in monthday)
        with xr.Dataset(data_vars=data_vars, attrs=attrs) as obj:
            obj.to_netcdf(prefix_out + datestr + '.nc', encoding=encoding)


########
# Perform the extraction
########

# Define indices and variables
indices = {'gridX': slice(0, 398), 'gridY': slice(0, 898), 'depth': 0}
variables = ['uVelocity', 'vVelocity']

dateranges = [
    (datetime(2014, 1, 1), datetime(2014, 3, 16)),
    (datetime(2014, 5, 1), datetime(2014, 7, 14)),
    (datetime(2014, 8, 1), datetime(2014, 10, 14)),
]

for daterange in dateranges:
   ds= extract_hindcast(daterange, variables, indices=indices)
