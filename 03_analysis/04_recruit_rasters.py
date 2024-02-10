# create raster of the count of points that strand or are over the 30m bathymetry
# at the end of a PLD
# This will be recruitment.


import arcpy
from arcpy.sa import *
import os


root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts'
destpts_dir = 'sim{}/outputs/shp/dest_pts'
output_gdb = 'DEST_RAST.gdb'
plds = [1,3,7,10,21,30,40,60]
sims = [
    '1101','1105','1108',
    '1401','1405','1408',
    '1701','1705','1708'
]

# get paths of dest_pts
arcpy.env.workspace = os.path.join(root, 'COMBINED.gdb')
files = arcpy.ListFeatureClasses(wild_card='destpts_*')
fcs = []
for fc in files:
    fcs.append(os.path.join(root, 'COMBINED.gdb', fc))

# create raster output gdb and set as workspace
if not arcpy.Exists(os.path.join(root, output_gdb)):
    arcpy.CreateFileGDB_management(root, output_gdb)
arcpy.env.workspace = os.path.join(root, output_gdb)


# copy to outgdb and Delete Identical points
# Where MPAs overlap, points get duplicated, which is what I needed to create the
# correct connection lines, but I don't want to double count them here.
for fc in fcs:
    outfeat = 'delidentical_'+ os.path.basename(fc)
    print ('deleting identical points for {}'.format(os.path.basename(fc)))
    arcpy.CopyFeatures_management(fc, outfeat)
    arcpy.DeleteIdentical_management(outfeat, ['traj_id', 'uID_part'])


# The spatial selections take way too long with 3 million points
# Originally I would do the spatial selects first so that I can do different attribute selections
# based on what people might want to see in the future, but for now I think I should
# just do what is quick. It is unlikely that anyone will want to see self recruitment
# in the raster outputs.
fcs = arcpy.ListFeatureClasses(wild_card='delidentical*')
for fc in fcs:
    print('selecting by attribute for {}'.format(fc))
    # pld = int(fc.split('_')[3][3:])
    # pld_int = int((pld * 24) / time_step_output)
    arcpy.MakeFeatureLayer_management(fc, 'temp_lyr')
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'NEW_SELECTION', '"uID_part" <> "dest_id"')
    # so the way I have dest_pts now, the max that time_int can be is the pld_int. This would mean that the particle was still drifting at the end of that PLD.
    # so now remove for mortality, but only ones that have their mortality step before they stranded. I'm allowing ones to live if they strand  on the coast before mortality selection. However, for some species, if stranding is too early, they still may need time to develop, so in the future I could try a scenario where mortality could still apply.
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"mortstep" < "time_int" and "mortstep" <> -1')
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"uID_part" IN (169, 170)') # the particles for these features are all messed up. Not sure what happened. Best to just remove it.
    arcpy.CopyFeatures_management('temp_lyr', 'selattr_'+ fc.split('_', 1)[1])
    arcpy.Delete_management('temp_lyr')


# features to use for spatial select
landmask = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_FINAL'
mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'
contours = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\B07_bathy_merge'
con30 = arcpy.SelectLayerByAttribute_management(contours, 'NEW_SELECTION', '"ContourMax" = 30')
if not arcpy.Exists('landmask_buff30'):
    arcpy.Buffer_analysis(landmask, 'landmask_buff30', 30)
if not arcpy.Exists('land_mpas_con30'):
    arcpy.Merge_management(['landmask_buff30', mpas, con30], 'land_mpas_con30')
arcpy.Delete_management('landmask_buff30')


# I've now changed the spatial selection to be a just a clip
fcs = arcpy.ListFeatureClasses(wild_card='selattr_*')
# select points by location
for fc in fcs:
    print('clipping to land_mpas_con30 for {}'.format(fc))
    # old way (hold on to this for a while):
    # arcpy.MakeFeatureLayer_management(fc, 'temp_lyr')
    # arcpy.SelectLayerByLocation_management('temp_lyr', 'INTERSECT', landmask, search_distance=50, selection_type='NEW_SELECTION')
    # arcpy.SelectLayerByLocation_management('temp_lyr', 'INTERSECT', mpas, selection_type='ADD_TO_SELECTION')
    # arcpy.SelectLayerByLocation_management('temp_lyr', 'INTERSECT', 'contour30', selection_type='ADD_TO_SELECTION')
    # arcpy.CopyFeatures_management('temp_lyr', 'selspatial_'+ fc.split('_', 1)[1])
    # arcpy.Delete_management('temp_lyr')
    arcpy.Clip_analysis(fc, 'land_mpas_con30', 'clip_'+ fc.split('_', 1)[1])



# so for a pld of 60, once I select for ones that didn't die, this pretty much removes everything
# that didn't strand, which means the spatial select by mpas and contour are not that meaningful.
# For a PLD of 60-75, with 1.6million particles, only 8 particles would be remaining at the end (with mortality).
# This sounds extreme, but it makes sense. Therefore, settlement is super important.
# I will see different patterns with shorter PLDs. Some in open water may still be living.
# This is also a good reason to perhaps have a mortality rate that changes through time,
# then I am not losing so many right away.


# convert to raster
snapras = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_NN'
arcpy.env.snapRaster = snapras
fcs = arcpy.ListFeatureClasses(wild_card='clip*')
# snap to the landmask raster which is 500m cell size. However, only snap to it the first time
# and then change the snap raster to the first one created. I am creating rasters with 1000m
# cell size, so if I keep snapping to the landmask then they don't always get created on
# the same grid
newSnapRas = True
for fc in fcs:
    pld = fc.split('_')[-1]
    date = fc.split('_')[-2]
    outname = 'recruit_count_{}_{}'.format(date, str(pld))
    arcpy.PointToRaster_conversion(fc, '', outname, 'COUNT', cellsize=1000)
    if newSnapRas:
        arcpy.env.snapRaster = outname
        newSnapRas = False


# add rasters by PLD to create one overall raster by PLD
rasters = arcpy.ListRasters(wild_card='recruit_count*')
for pld in plds:
    ras_plds = []
    for ras in rasters:
        base_pld = ras.split('_')[-1]
        if base_pld == 'pld{}'.format(str(pld)):
            ras_plds.append(ras)
    outras = CellStatistics(ras_plds, 'SUM', 'DATA')
    outras.save('recruit_count_ALL_{}'.format(str(pld)))


#
# remove cells that overlap the BUFFERED MPA (turn to null values)
# Then I can symbolize the unprotected areas better.
#

# rasterize mpa buff
mpa_buff = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpas.gdb\MB06_mpa_buff_FINAL'
mpa_buff_out = os.path.basename(mpa_buff)
arcpy.CopyFeatures_management(mpa_buff, mpa_buff_out)
arcpy.AddField_management(mpa_buff_out, 'priority', 'SHORT')
with arcpy.da.UpdateCursor(mpa_buff_out, ['priority']) as cursor:
    for row in cursor:
        row[0]=1
        cursor.updateRow(row)
snapras = 'recruit_count_ALL_{}'.format(str(plds[0]))
arcpy.env.snapRaster = snapras
arcpy.env.cellSize = snapras
arcpy.PolygonToRaster_conversion(mpa_buff_out, 'OBJECTID', mpa_buff_out + '_ras', 'MAXIMUM_AREA', priority_field='priority')
# its super important to have a priority field and the conversion set to maximum area or else it won't give a cell a value unless the polygon crosses the cell center

# give the cells a value of 0 where there is an mpa and a value of 1 elsewhere:
outCon = Con(mpa_buff_out + '_ras', 0, 1, "Value > 0")
outisnull = IsNull(outCon)
outisnull.save(mpa_buff_out + '_ras_isnull')

# use the set null tool
# If mparast is zero then the output raster will be null there, but where it is not then it will use the values from recruitras
rasters = arcpy.ListRasters(wild_card='recruit_count_ALL_*')
mpas_buff_isnull = mpa_buff_out + '_ras_isnull'
for ras in rasters:
    outras = SetNull(mpas_buff_isnull, ras, 'Value = 0')
    outras.save(ras + '_rmMPAs')


