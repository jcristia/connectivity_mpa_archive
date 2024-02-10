# prepare coastline to act as landmask in OpenDrift

import arcpy
import netCDF4 as nc
import numpy as np

# inputs
coast_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb'
arcpy.env.workspace = coast_gdb

nep36 = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\mesh_mask.nc"
out_csv = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\coastline_mpas\mask_nepnemo_JC.csv'

dataset = nc.Dataset(nep36, "r+")


# # explore
# print(dataset.data_model)
# variables = dataset.variables.keys()
# for variable in variables:
#     print(variable)
# for dim in dataset.dimensions.values():
#     print(dim)
# for var in dataset.variables.values():
#     print(var)
# tmask = dataset.variables["tmask"][:]
# np.shape(tmask)
# tmask[0][0] 
# tmask[:,0]


# extract xy and landmask
lat = dataset.variables["nav_lat"][:]
lon = dataset.variables["nav_lon"][:]
tmask = dataset.variables['tmask'][:]
tmask = tmask[0][0]
lat = np.array(lat).flatten()
lon = np.array(lon).flatten()
tmask = np.array(tmask).flatten()
# switch land and water
tmask = 1-tmask

headers = ['X', 'Y', 'land']
output = open(out_csv, 'w')
for header in headers:
    output.write(header + ",")
output.write("\n")
for x,y,tmask in zip(lon,lat, tmask):
    if x != 0 and y != 0:
        # can only write strings
        output.write('%s' % x)
        output.write(',')
        output.write('%s' % y)
        output.write(',')
        output.write('%s' % tmask)
        output.write(',')
        output.write("\n")
output.close()


# create point dataset
arcpy.management.XYTableToPoint(out_csv, "landmask_pt_wgs84", 'X', 'Y', coordinate_system=arcpy.SpatialReference(4326))
# project to BC albers
arcpy.Project_management('landmask_pt_wgs84', 'landmask_pt_BCalbers', arcpy.SpatialReference(3005))

# nearest neighbor interpolation
# NN interpolation does the best job of maintaining the structure of the land/water point structure. "Its basic properties are that it's local, using only a subset of samples that surround a query point, and interpolated heights are guaranteed to be within the range of the samples used. It does not infer trends and will not produce peaks, pits, ridges, or valleys that are not already represented by the input samples."
arcpy.NaturalNeighbor_3d('landmask_pt_BCalbers', 'land', 'landmask_NN', 500)


# contour
arcpy.Contour_3d('landmask_NN', 'landmask_contour', 0.1, 0, 1, 'CONTOUR_SHELL_UP')

# extract just the level I need
# Based on visual assessment, select 0.7 level
# select by attribute and save to new fc
lm = arcpy.SelectLayerByAttribute_management('landmask_contour', 'NEW_SELECTION', '"ContourMin" = 0.7')
arcpy.CopyFeatures_management(lm, 'landmask_contour_07')






###################################################
# NOW, do this all for the SalishSeaCast
ssc_mask = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\salishsea\ubcSSn3DMeshMaskV17-02.nc"
ssc_latlon = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\salishsea\ubcSSnBathymetryV17-02.nc"
out_csv = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\coastline_mpas\mask_ssc_JC.csv'

dataset_m = nc.Dataset(ssc_mask, "r+")
dataset_l = nc.Dataset(ssc_latlon, "r+")


# extract xy and landmask
lat = dataset_l.variables["latitude"][:]
lon = dataset_l.variables["longitude"][:]
tmask = dataset_m.variables['tmask'][:]
tmask = tmask[0][0]
lat = np.array(lat).flatten()
lon = np.array(lon).flatten()
tmask = np.array(tmask).flatten()
# switch land and water
tmask = 1-tmask

headers = ['X', 'Y', 'land']
output = open(out_csv, 'w')
for header in headers:
    output.write(header + ",")
output.write("\n")
for x,y,tmask in zip(lon,lat, tmask):
    if x != 0 and y != 0:
        # can only write strings
        output.write('%s' % x)
        output.write(',')
        output.write('%s' % y)
        output.write(',')
        output.write('%s' % tmask)
        output.write(',')
        output.write("\n")
output.close()

arcpy.management.XYTableToPoint(out_csv, "landmask_pt_wgs84_ssc", 'X', 'Y', coordinate_system=arcpy.SpatialReference(4326))
arcpy.Project_management('landmask_pt_wgs84_ssc', 'landmask_pt_BCalbers_ssc', arcpy.SpatialReference(3005))
arcpy.NaturalNeighbor_3d('landmask_pt_BCalbers_ssc', 'land', 'landmask_NN_ssc', 100)
arcpy.Contour_3d('landmask_NN_ssc', 'landmask_contour_ssc', 0.1, 0, 1, 'CONTOUR_SHELL_UP')
lm = arcpy.SelectLayerByAttribute_management('landmask_contour_ssc', 'NEW_SELECTION', '"ContourMin" = 0.7')
arcpy.CopyFeatures_management(lm, 'landmask_contour_ssc_07')



################################################

# combine the two
# get minimum bounding rectangle to erase with
arcpy.MinimumBoundingGeometry_management('landmask_contour_ssc_07', 'landmask_mb_ssc', 'CONVEX_HULL')
# !!!!!! MANUAL EDIT:
# From this bounding rectangle AND the ssc landmask, remove barkley sound and alberni inlet
arcpy.CopyFeatures_management('landmask_mb_ssc', 'landmask_mb_ssc_manualedit')
arcpy.CopyFeatures_management('landmask_contour_ssc_07', 'landmask_contour_ssc_07_manualedit')

arcpy.Erase_analysis('landmask_contour_07', 'landmask_mb_ssc_manualedit', 'landmask_contour_07_erase')
arcpy.Merge_management(['landmask_contour_07_erase','landmask_contour_ssc_07_manualedit'], 'landmask_merge')
arcpy.Dissolve_management('landmask_merge','landmask_dissolve', multi_part='SINGLE_PART')


# !!!!! MANUAL EDITING HAPPENING HERE
# I need to remove some slivers, which I don't want to automate
# and I need to delete the open boundary lines for the ssc data
# editing will happen on this copy:
arcpy.CopyFeatures_management('landmask_dissolve', 'landmask_manualedit')

# convert to shapefile
arcpy.CopyFeatures_management('landmask_manualedit', 'landmask_FINAL')
arcpy.FeatureClassToShapefile_conversion('landmask_FINAL', r"C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline")

# project (opendrift needs WGS84, but I still want the BC Albers version for mapping)
arcpy.Project_management(r"C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\landmask_FINAL.shp", r"C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\landmask_FINAL_wgs84.shp", arcpy.SpatialReference(4326))