# average connectivity over all time periods
# (but not over all PLDs. I won't do that for this chapter)


import os
import pandas as pd
import geopandas as gp

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts'
shps = os.path.join(root, 'sim{}/outputs/shp')
sims = [
    '1101','1105','1108',
    '1401','1405','1408',
    '1701','1705','1708'
]
plds = [1,3,7,10,21,30,40,60]


def mean_cust_denom(x):
    s = x.sum()
    m = s/float(len(sims))
    return m


for pld in plds:

    # get all the shapefiles for a pld across all times
    files = []
    for sim in sims:
        folder = shps.format(sim)
        folder = os.path.join(folder, 'conn_lines')
        fl = os.listdir(folder)
        for file in fl:
            base_pld = file.split('_')[-1]
            if base_pld == 'pld{}.shp'.format(str(pld)):
                files.append(os.path.join(folder, file))

    # put them all in geodataframe
    gdf_all = gp.GeoDataFrame()
    for shp in files:
        gdf = gp.read_file(shp)
        gdf_all = gdf_all.append(gdf)
    gdf_all['pld'] = int(pld)
    gdf_all = gdf_all.astype({'from_id':int, 'to_id':int}) # there's still a mix of datatypes in the columns for some reason. This was super important to do or else the code below didn't recognize duplicates.

    # average
    gdf_group = gdf_all.groupby(['from_id', 'to_id']).agg(
        freq = ('from_id', 'count'),
        prob_avg = ('prob', mean_cust_denom),
        #time_int = ('time_int', 'first'), # no longer including this since it was not included correctly from the biology script
        totalori = ('totalori', 'sum'),
        totquant = ('quantity', 'sum'),
        date_start = ('date_start', 'first'),
        geometry = ('geometry', 'first'),
        pld = ('pld', 'first')
        )
    gdf_group = gdf_group.astype({'totalori':int, 'totquant':int})
    gdf_group = gdf_group.reset_index()

    # output
    gdf_f = gp.GeoDataFrame(gdf_group, crs=gdf.crs)
    gdf_f.to_file(filename=os.path.join(root, 'conn_avg_pld{}.shp'.format(str(pld))), driver='ESRI Shapefile')


# output to feature class
# need to switch to arcpro environment
# kinda hokey, but you can't install geopandas in the arcpro environment
# delete shapefiles after
import arcpy
root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts'
gdb = r'COMBINED.gdb'
arcpy.env.workspace = os.path.join(root, gdb)
fl = os.listdir(root)
for file in fl:
    name = os.path.splitext(file)[0]
    ext = os.path.splitext(file)[1]
    if ext == '.shp':
        arcpy.CopyFeatures_management(os.path.join(root, file), name)

# delete shapefiles
fl = os.listdir(root)
for file in fl:
    name = os.path.splitext(file)[0]
    ext = os.path.splitext(file)[1]
    if ext == '.shp':
        arcpy.Delete_management(os.path.join(root, file))