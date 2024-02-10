# create a buffered version of the mpas dataset
# this is used for settlement in the biology module
# particles that strand may be partially inland
# a spot check of a few destination points datasets shows that a 100m
# buffer should be adequate

import arcpy

mpa_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_buff.gdb'
arcpy.env.workspace = mpa_gdb

mpas_og = r'M08_mpa_20201124_FINALCOPY'
landmask = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_FINAL'

mpas_shp_folder_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpas.gdb'

# select by location those that touch the coastline
# (I don't need to buffer all the small ones that don't touch land)
sel_intersect = arcpy.SelectLayerByLocation_management(mpas_og, 'INTERSECT', landmask, selection_type='NEW_SELECTION')
sel_fc = arcpy.CopyFeatures_management(sel_intersect, 'MB01_mpa_selectbylocation')

# buffer by 100m
buff_fc = arcpy.Buffer_analysis(sel_fc, 'MB02_mpa_buffer100m', '200 meters', dissolve_option='NONE')

# clip to coastline
clipped_fc = arcpy.Clip_analysis(buff_fc, landmask, 'MB03_mpa_cliptoland')

# merge this with original
merged_fc = arcpy.Merge_management([mpas_og, clipped_fc], 'MB04_mpa_merged')

# dissolve by uID_part
# (this may add donuts, but I don't think this matters for the biology script)
# Also overlap of adjacent MPAs is ok now.
diss_fc = arcpy.Dissolve_management(merged_fc, 'MB05_mpa_dissolve', 'uID_20201124')

# copy final buff fc to mpas_shp_release gdb
final_copy = arcpy.CopyFeatures_management(diss_fc, 'MB06_mpa_buff_FINAL')
arcpy.CopyFeatures_management(final_copy, os.path.join(mpas_shp_folder_gdb, 'MB06_mpa_buff_FINAL'))

# then manually copy it to shp outside of gdb
# C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\