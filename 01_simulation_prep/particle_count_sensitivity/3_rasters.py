# Convert points in rasters.
# Following Simons et al. 2013...
# create rasters based on a random selection of x% of points. 5-100%
# Compare all the rasters to the raster that includes all the points. Use the
# band statistics tool in ArcGIS to obatin a correlation matrix.


mpa_uid = 207


import arcpy
from arcpy.sa import *
import os
import random

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particle_count\analysis'
shp_dir = os.path.join(root, 'dest_pt_shp')
rst_dir = os.path.join(root, 'dest_rast')
bnd_dir = os.path.join(root, 'band_statistics')

gdb = f'dest_rast_{str(mpa_uid)}.gdb'
arcpy.env.workspace = os.path.join(rst_dir, gdb)

# create gdb for each mpa
if not arcpy.Exists(arcpy.env.workspace):
    arcpy.CreateFileGDB_management(rst_dir, gdb)

# create sim number attribute
for file in os.listdir(shp_dir):
    fn, fext = os.path.splitext(file)
    mid = int(fn.split('_')[1])
    sim = int((fn.split('_')[2]).split('.')[0])
    if fext=='.shp' and mid==mpa_uid:
        shp_path = os.path.join(shp_dir, file)        
        arcpy.AddField_management(shp_path, 'simulation', 'SHORT')
        arcpy.CalculateField_management(shp_path, 'simulation', sim)

# add shapefile to list and merge all
shp_pts = []
for file in os.listdir(shp_dir):
    fn, fext = os.path.splitext(file)
    mid = int(fn.split('_')[1])
    sim = int((fn.split('_')[2]).split('.')[0])
    if fext=='.shp' and mid==mpa_uid:
        shp_path = os.path.join(shp_dir, file)        
        shp_pts.append(shp_path)
arcpy.Merge_management(shp_pts, 'dest_pts_ALL')



#### create 100% raster
# raster requires value field. Add dummy field with value of 1
arcpy.AddField_management('dest_pts_ALL', 'value', 'SHORT')
arcpy.CalculateField_management('dest_pts_ALL', 'value', 1)
# create raster
arcpy.PointToRaster_conversion('dest_pts_ALL', 'value', 'dest_rast_100TEMP', 'COUNT', cellsize=5000)
# divide by total points
pt_count = int(arcpy.GetCount_management('dest_pts_ALL')[0])
ras = Raster('dest_rast_100TEMP')
floatras = ras / pt_count
outIsNull = IsNull(floatras)
outCon = Con(outIsNull, 0, floatras)
outCon.save('dest_rast_100')
arcpy.Delete_management('dest_rast_100TEMP')

# NOTE: it was very important to turn nodata cells to 0. I was getting values
# less than -1 in the correlation matrix, which should not be possible.


#### create rasters from random selection of points
# do for 20 levels
arcpy.env.snapRaster = 'dest_rast_100'
arcpy.env.cellsize = 'dest_rast_100'
arcpy.env.extent = 'dest_rast_100'
for i in range(5, 100, 5):
    pt_sample = round((i/100.0) * pt_count)
    rnd_samp = random.sample(range(1,pt_count), pt_sample)
    
    # select from dest_pts_ALL and make new raster
    sel = arcpy.SelectLayerByAttribute_management(
        'dest_pts_ALL',
        'NEW_SELECTION',
        '"traj_id" IN {}'.format(str(tuple(rnd_samp)))
    )
    arcpy.PointToRaster_conversion(sel, 'value', f'dest_rast_{i}TEMP', 'COUNT')

    # divide by total points
    ras = Raster(f'dest_rast_{i}TEMP')
    floatras = ras / pt_count
    outIsNull = IsNull(floatras)
    outCon = Con(outIsNull, 0, floatras)
    outCon.save(f'dest_rast_{i}')
    arcpy.Delete_management(f'dest_rast_{i}TEMP')


# compare each level to overall raster
rasters = []
for i in range(5, 105, 5):
    r = f'dest_rast_{i}'
    rasters.append(r)
outStatFile = os.path.join(bnd_dir, f'bandstats_mpa_{mpa_uid}.txt')
BandCollectionStats(rasters, outStatFile, 'DETAILED')
