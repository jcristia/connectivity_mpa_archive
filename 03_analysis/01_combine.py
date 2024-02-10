# combine shapefile outputs for biology script

# for each time period, and for each of 9 groups:
# there will be connection line shapefiles for each of 8 PLDs
# there will be 1 destination point shapefiles for the 60 day PLD


import arcpy
import os


root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts'
shps = os.path.join(root, 'sim{}/outputs/shp')
sims = [
    '1101','1105','1108',
    '1401','1405','1408',
    '1701','1705','1708'
]
plds = [1,3,7,10,21,30,40,60]
gdb = r'COMBINED.gdb'
arcpy.env.workspace = os.path.join(root, gdb)
# create gdb if it doesn't exist
if not arcpy.Exists(arcpy.env.workspace):
    arcpy.CreateFileGDB_management(root, gdb)


# copy in centroids once
if not arcpy.Exists('patch_centroids'):
    centroids = os.path.join(shps.format(sims[0]), 'patch_centroids.shp')
    arcpy.CopyFeatures_management(centroids, 'patch_centroids')


# for each time period
# combine the line shapefiles for each group
# output to fgdb so that they can all be in 1 place
for sim in sims:
    folder = shps.format(sim)
    folder = os.path.join(folder, 'conn_lines')
    fl = os.listdir(folder)
    for pld in plds:
        files = []
        for file in fl:
            base_pld = file.split('_')[-1]
            if base_pld == 'pld{}.shp'.format(str(pld)):
                files.append(os.path.join(folder, file))
        fc = 'connectivity_{}_pld{}'.format(sim, str(pld))
        if not arcpy.Exists(fc):
            arcpy.Merge_management(files, fc)

# for some reason from_id is a text field
for fc in arcpy.ListFeatureClasses('connectivity*'):
    arcpy.AddField_management(fc, 'f_temp', 'SHORT')
    #arcpy.CalculateField_management(fc, 'f_temp', '!from_id!')
    # Calculate field is giving me a weird error that Python is not installed.
    with arcpy.da.UpdateCursor(fc, ['from_id', 'f_temp']) as cursor:
        for row in cursor:
            row[1] = int(row[0])
            cursor.updateRow(row)
    arcpy.DeleteField_management(fc, ['from_id'])
    arcpy.AddField_management(fc, 'from_id', 'SHORT')
    #arcpy.CalculateField_management(fc, 'from_id', '!f_temp!')
    with arcpy.da.UpdateCursor(fc, ['f_temp', 'from_id']) as cursor:
        for row in cursor:
            row[1] = row[0]
            cursor.updateRow(row)
    arcpy.DeleteField_management(fc, ['f_temp'])


# merge destination points per time period
# this is going to be 3 million points per fc
for sim in sims:
    folder = shps.format(sim)
    folder = os.path.join(folder, 'dest_pts')
    fl = os.listdir(folder)
    for pld in plds:
        files = []
        for file in fl:
            base_pld = file.split('_')[-1]
            if base_pld == 'pld{}.shp'.format(str(pld)):
                files.append(os.path.join(folder, file))
        fc = 'destpts_{}_pld{}'.format(sim, str(pld))
        if len(files) > 0 and not arcpy.Exists(fc):
            arcpy.Merge_management(files, fc)
