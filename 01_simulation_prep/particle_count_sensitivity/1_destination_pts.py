
# For each mpa's netcdf files, get the destination pts into feature classes.


mpa_id = 207
sims = 20

import netCDF4 as nc
import os
import numpy as np
import pandas as pd
from shapely.geometry import shape, Point, LineString, Polygon
import geopandas as gp
import logging
logging.basicConfig(level=logging.INFO)


nc_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particle_count\simulations\outputs\nc'
out_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particle_count\analysis\dest_pt_shp'


def get_destination_coords(traj, lon, lat, status):

    logging.info("getting destination coordinates")
    lons_dest = []
    lats_dest = []
    for i in range(len(traj)):
        if np.ma.is_masked(lon[i][-1]): # if the last value is masked (but really just to check if any values are masked, similar to getParticleOriginPoly). If it is masked then it must have stranded, and therefore we can search by where it is 1.
            # changed this statement from '== 1' to '> 0'. I found that if a particle goes outside of the grid it gets coded as '2 - missing data'. Opendrift says that anything above 0 is considered deactivated, so the actual number doesn't matter (at least for my purposes, just that it is bigger than 0.
            j = np.where(status[i] > 0)[0][0]
            lo = lon[i][j]
            lons_dest.append(lo)
            la = lat[i][j]
            lats_dest.append(la)
        else: # otherwise just get the last coordinate
            lons_dest.append(lon[i][-1])
            lats_dest.append(lat[i][-1])

    df = pd.DataFrame()
    df['Coordinates'] = list(zip(lons_dest, lats_dest))
    df['Coordinates'] = df['Coordinates'].apply(Point)
    df['traj_id'] = list(traj)
    points_dest = gp.GeoDataFrame(df, geometry='Coordinates')
    points_dest.crs = {'init' :'epsg:4326'}
    points_dest = points_dest.to_crs(epsg=3005)
    points_dest = points_dest.infer_objects()
    points_dest.traj_id = points_dest.traj_id.astype('float')

    return points_dest




def out_shp_dest_points(points_dest, shp_out):

    points_dest.to_file(filename=shp_out, driver='ESRI Shapefile')


ncs = []
for file in os.listdir(nc_dir):
    mid = int(file.split('_')[1])
    sim = int((file.split('_')[2]).split('.')[0])
    if mid == mpa_id and sim < sims:
        ncs.append(file)

for n in ncs:

    mid = int(n.split('_')[1])
    sim = int((n.split('_')[2]).split('.')[0])
    shp_out = os.path.join(out_dir, f'destpts_{mid}_{sim}.shp')

    dataset = nc.Dataset(os.path.join(nc_dir, n), "r+")
    lon = dataset.variables["lon"]
    lat = dataset.variables["lat"]
    traj = dataset.variables["trajectory"]
    status = dataset.variables["status"]

    points_dest = get_destination_coords(traj, lon, lat, status)
    out_shp_dest_points(points_dest, shp_out)