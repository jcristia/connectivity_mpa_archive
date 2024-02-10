# prep the MPA dataset
# create bathymetry contours
# some manual editing required

import arcpy
import netCDF4 as nc
import numpy as np

mpa_gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb'
arcpy.env.workspace = mpa_gdb

# Canada and USA mpa datasets
mpa_can = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_protected_areas.gdb\CPCAD_Dec2019'
mpa_usa = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_protected_areas.gdb\NOAA_MPAI_2020_IUCN'

# hydrodynamic models to create bathymetry contours
nep36 = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\nep_nemo\Bathymetry_NEP36_714x1020_SRTM30v11_NOAA3sec_WCTSS_JdeFSalSea.nc"
ssc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\spatial\models\salishsea\ubcSSnBathymetryV17-02.nc'

# coastlines used for clipping
# these need to be detailed, but it doesn't matter too much since I will need to fill in holes and align the mpas to the generalized landmask anyways
coast_canada = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_boundaries.gdb\Provinces'
coast_usa = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_boundaries.gdb\states_detailed'
nep36_mesh = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_pt_BCalbers'  # the points for creating a minimum bounding area for clipping
ssc_mask = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\Coastline\coastline.gdb\landmask_mb_ssc_manualedit'



##############################################


# create a minimum bounding polygons
arcpy.MinimumBoundingGeometry_management(nep36_mesh, 'A00_nep36_mb', 'CONVEX_HULL')
arcpy.CopyFeatures_management(ssc_mask, 'A00_ssc_mb') # used later for contours

# clip each mpa dataset
arcpy.Clip_analysis(mpa_can, 'A00_nep36_mb', 'M01_mpa_can_clip')
arcpy.Clip_analysis(mpa_usa, 'A00_nep36_mb', 'M01_mpa_usa_clip')

# erase mpas with coastlines
arcpy.Erase_analysis('M01_mpa_can_clip', coast_canada, 'M02_mpa_can_erase')
arcpy.Erase_analysis('M01_mpa_usa_clip', coast_usa, 'M02_mpa_usa_erase')

# add field for country
arcpy.AddField_management('M02_mpa_can_erase', 'country', 'TEXT')
arcpy.AddField_management('M02_mpa_usa_erase', 'country', 'TEXT')
arcpy.CalculateField_management('M02_mpa_can_erase', 'country', '"Canada"', 'PYTHON_3')
arcpy.CalculateField_management('M02_mpa_usa_erase', 'country', '"USA"', 'PYTHON_3')

# Select and copy Marine features for CARTS
marine = arcpy.SelectLayerByAttribute_management('M02_mpa_can_erase', 'NEW_SELECTION', "BIOME = 'M'") # must use Biome. Other descriptions are not complete
arcpy.CopyFeatures_management(marine, 'M03_mpa_can_selmarine')

# merge
arcpy.Merge_management(['M03_mpa_can_selmarine', 'M02_mpa_usa_erase'], 'M04_mpa_merge')

# create a copy
arcpy.CopyFeatures_management('M04_mpa_merge', 'M05_mpa_MASTER')
# add unique ID and calculate
arcpy.AddField_management('M05_mpa_MASTER', 'uID', 'SHORT')
arcpy.CalculateField_management('M05_mpa_MASTER', 'uID', '!OBJECTID!', 'PYTHON_3')
# (this will be the master version to join back to)



#################################################
# create bathymetry contours
ds_nep36 = nc.Dataset(nep36, "r+")
ds_ssc = nc.Dataset(ssc, "r+")
out_csv_nep = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\coastline_mpas\bathy_nepnemo_JC.csv'
out_csv_ssc = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\coastline_mpas\bathy_ssc_JC.csv'
# extract xy and bathy for nep 36
lat = ds_nep36.variables["nav_lat"][:]
lon = ds_nep36.variables["nav_lon"][:]
bath = ds_nep36.variables['Bathymetry'][:]
lat = np.array(lat).flatten()
lon = np.array(lon).flatten()
bath = np.array(bath).flatten()
headers = ['X', 'Y', 'bathy']
output = open(out_csv_nep, 'w')
for header in headers:
    output.write(header + ",")
output.write("\n")
for x,y,bath in zip(lon,lat, bath):
    if x != 0 and y != 0:
        # can only write strings
        output.write('%s' % x)
        output.write(',')
        output.write('%s' % y)
        output.write(',')
        output.write('%s' % bath)
        output.write(',')
        output.write("\n")
output.close()

# extract xy and bathy for ssc
lat = ds_ssc.variables["latitude"][:]
lon = ds_ssc.variables["longitude"][:]
bath = ds_ssc.variables['bathymetry'][:]
lat = np.array(lat).flatten()
lon = np.array(lon).flatten()
bath = np.array(bath).flatten()
headers = ['X', 'Y', 'bathy']
output = open(out_csv_ssc, 'w')
for header in headers:
    output.write(header + ",")
output.write("\n")
for x,y,bath in zip(lon,lat, bath):
    if x != 0 and y != 0:
        # can only write strings
        output.write('%s' % x)
        output.write(',')
        output.write('%s' % y)
        output.write(',')
        output.write('%s' % bath)
        output.write(',')
        output.write("\n")
output.close()

# create point dataset
arcpy.management.XYTableToPoint(out_csv_nep, "B01_bathy_pt_wgs84_nep36", 'X', 'Y', coordinate_system=arcpy.SpatialReference(4326))
arcpy.management.XYTableToPoint(out_csv_ssc, "B01_bathy_pt_wgs84_ssc", 'X', 'Y', coordinate_system=arcpy.SpatialReference(4326))
# project to BC albers
arcpy.Project_management('B01_bathy_pt_wgs84_nep36', 'B02_bathy_pt_BCalbers_nep36', arcpy.SpatialReference(3005))
arcpy.Project_management('B01_bathy_pt_wgs84_ssc', 'B02_bathy_pt_BCalbers_ssc', arcpy.SpatialReference(3005))
# nearest neighbor interpolation
arcpy.NaturalNeighbor_3d('B02_bathy_pt_BCalbers_nep36', 'bathy', 'B03_bathy_NN_nep', 500)
arcpy.NaturalNeighbor_3d('B02_bathy_pt_BCalbers_ssc', 'bathy', 'B03_bathy_NN_ssc', 100)
# contour
arcpy.Contour_3d('B03_bathy_NN_nep', 'B04_bathy_contour_nep', 15, 0, 1, 'CONTOUR_SHELL')
arcpy.Contour_3d('B03_bathy_NN_ssc', 'B04_bathy_contour_ssc', 15, 0, 1, 'CONTOUR_SHELL')

# my justification for using 30m:
# yes, we are doing surface dispersal, but given the resolution of the models, we need to generalize and create strips of MPAs large enough so that we can give particles the chance to release into open water. 30m is just beyond the max depth of seagrass and kelp, so we aren't generalizing too much

# export out 15-45m bathymetry
sel_nep = arcpy.SelectLayerByAttribute_management('B04_bathy_contour_nep', 'NEW_SELECTION', "ContourMax <= 45")
arcpy.CopyFeatures_management(sel_nep, 'B05_bathy_1545_nep')
sel_ssc = arcpy.SelectLayerByAttribute_management('B04_bathy_contour_ssc', 'NEW_SELECTION', "ContourMax <= 45")
arcpy.CopyFeatures_management(sel_ssc, 'B05_bathy_1545_ssc')
# erase salish sea area from nep
arcpy.Erase_analysis('B05_bathy_1545_nep', 'A00_ssc_mb', 'B06_bathy_erase_nep')
# merge
arcpy.Merge_management(['B05_bathy_1545_ssc', 'B06_bathy_erase_nep'], 'B07_bathy_merge')



########################################

# MANUAL WORK

# Instead of editing the exisiting mpa dataset, I will just draw new polygons. This will actually be simpler since my land polys are so simple.

# Manually create a new feature class with todays date.
# Attributes: uID_{DATE}, orig_uID, partID, comments  

#### Feature creation instructions: ####
# in the master mpa datset, go line by line
# have open for reference: my mpa dataset, bathymetry (0-45m) contours, landmask coastline, original coastline, pt bathymetry values, Katie's MPA dataset
# when viewing an mpa, be aware of multipart features. Make sure you are viewing the entire thing.
# Draw them as single part features, give them the same orig_uID but different partID numbers
# When drawing MPAs around entire islands, draw it so that it is technically not a hole (Opendrift requirements)
# Look for obviously nearshore MPAs OR MPAs that overlap with 30m bathymetry
# Trace it so that it algins with the coast even if it doesn't perfectly align in the detailed file. The exception would be if it is one that overlaps only with submerged areas and is not anywhere near the coast.
# If it extends beyond 30m then trace it out just to 30m bathymetry.
# If it does not extend out to 30m, then extend it to just its actual extent.
# For cases where the 30m contour is on or super close to land: extend it out to the next contour that is completely off land for that strip of MPA. This seems mostly to occur in the NEP model. It looks like there are some extreme cases where I will need to keep the original bathymetry on in the background use, for example, a 150m contour. This would be for areas were the detailed coast extends out farther than the simplified coast AND the bottom drops off instantly.

# As I worked through this, I found that there are some small mpas around islets where the land does not show up in the generalized coastlin, and in some cases there is not land at all but they are within 20m so they show up in Katie's MPA dataset.
# In this case, I included some of them. I will just need to say that I used a reference 10m bathymetry for some small islet areas.
# They would have to include a significant amount of seafloor above 20m for me to include it.
# Since I am only modeling surface dispersal I don't want to include little specks.

########################################

# The last feature is the San Juan Islands Wildlife Refuge
# However, it is just the not take areas, which are a bunch of islets.
# It think it is important to capture these details in case I choose not to inlcude the rest of the general reserve, which doesn't seem like it has a high level of protection.

# Run these steps after doing the above manual work:

# select just the San Juan Islets feature
sel_sji = arcpy.SelectLayerByAttribute_management('M01_mpa_usa_clip', 'NEW_SELECTION', "Site_ID = 'NWR200'")
arcpy.CopyFeatures_management(sel_sji, 'S01_mpa_sji')
# clean and generalize
arcpy.MultipartToSinglepart_management('S01_mpa_sji', 'S02_mpa_sji_multising')
arcpy.AggregatePolygons_cartography('S02_mpa_sji_multising', 'S03_mpa_sji_aggregate', 250)
# Generalize shape (we don't need this amount of detail and it looks like weird almost holes are created)
arcpy.SimplifyPolygon_cartography('S03_mpa_sji_aggregate', 'S04_mpa_sji_simplify', 'POINT_REMOVE', 20, '', '', 'NO_KEEP')
# copy
arcpy.CopyFeatures_management('S04_mpa_sji_simplify', 'S05_mpa_sji_manualedit')
# !!! Manual edits !!!!
# in this dataset, remove any weird geometries that Opendrift may interpret as holes
# also there are 4 large-ish islands in this dataset. Remove their area.

# merge with working dataset (append requires working with attributes)
arcpy.Merge_management(['M06_mpa_20201124','S05_mpa_sji_manualedit'], 'M07_mpa_20201124_mergeS05')
# !!!!!! Manual Work !!!!!!
# fill in attributes for the new features
# !!!! I noticed that my unique IDs were off. However, after the merge my ObjectID resets and I can just use them to recalculate.
# then make a copy so that I can label it FINAL
arcpy.CopyFeatures_management('M07_mpa_20201124_mergeS05', 'M08_mpa_20201124_FINALCOPY')

# Create an archived version somewhere else for guaranteed safe keeping


############################################


# ATTRIBUTES and JOINING

# first, using the complete MPA dataset, create a new dataset with relevant attributes merged by field mapping

# split out by country (I guess I should have done the field mapping to begin with. Oh well.)
sel_can1 = arcpy.SelectLayerByAttribute_management('M05_mpa_MASTER', 'NEW_SELECTION', "country = 'Canada'")
sel_usa1 = arcpy.SelectLayerByAttribute_management('M05_mpa_MASTER', 'NEW_SELECTION', "country = 'USA'")
arcpy.CopyFeatures_management(sel_can1, 'M05_mpa_MASTER_tempCanada')
arcpy.CopyFeatures_management(sel_usa1, 'M05_mpa_MASTER_tempUSA')

can = 'M05_mpa_MASTER_tempCanada'
usa = 'M05_mpa_MASTER_tempUSA'
out_file = 'M05_mpa_MASTER_fieldsmapped'

# Fields (new name, canada name, usa name)
# uID: uID, uID
# govID: PARENT_ID, Site_ID
# name: NAME_E, Site_Name
# provstate: LOC_E, State
# gov_lev: GOV_TYPE, Gov_Level
# mgmt: MGMT_E, Mgmt_Agen
# type: TYPE_E, Design
# IUCN_CAT: IUCN_CAT, IUCNcat

fms = arcpy.FieldMappings()

fm_uID = arcpy.FieldMap()
fm_uID.addInputField(can, 'uID')
fm_uID.addInputField(usa, 'uID')
uID_name = fm_uID.outputField
uID_name.name = 'uID'
fm_uID.outputField = uID_name
fms.addFieldMap(fm_uID)

fm_govID = arcpy.FieldMap()
fm_govID.addInputField(usa, 'Site_ID')  # if merging text and numeric, put the text field first
fm_govID.addInputField(can, 'PARENT_ID')
govID_name = fm_govID.outputField
govID_name.name = 'govID'
govID_name.aliasName = 'govID'
fm_govID.outputField = govID_name
fms.addFieldMap(fm_govID)

fm_name = arcpy.FieldMap()
fm_name.addInputField(can, 'NAME_E')
fm_name.addInputField(usa, 'Site_Name')
name_name = fm_name.outputField
name_name.name = 'name'
name_name.aliasName = 'name'
fm_name.outputField = name_name
fms.addFieldMap(fm_name)

fm_ps = arcpy.FieldMap()
fm_ps.addInputField(can, 'LOC_E')
fm_ps.addInputField(usa, 'State')
ps_name = fm_ps.outputField
ps_name.name = 'provstate'
ps_name.aliasName = 'provstate'
fm_ps.outputField = ps_name
fms.addFieldMap(fm_ps)

fm_gl = arcpy.FieldMap()
fm_gl.addInputField(can, 'GOV_TYPE')
fm_gl.addInputField(usa, 'Gov_Level')
gl_name = fm_gl.outputField
gl_name.name = 'gov_lev'
gl_name.aliasName = 'gov_lev'
fm_gl.outputField = gl_name
fms.addFieldMap(fm_gl)

fm_man = arcpy.FieldMap()
fm_man.addInputField(can, 'MGMT_E')
fm_man.addInputField(usa, 'Mgmt_Agen')
man_name = fm_man.outputField
man_name.name = 'mgmt'
man_name.aliasName = 'mgmt'
fm_man.outputField = man_name
fms.addFieldMap(fm_man)

fm_type = arcpy.FieldMap()
fm_type.addInputField(can, 'TYPE_E')
fm_type.addInputField(usa, 'Design')
type_name = fm_type.outputField
type_name.name = 'type'
type_name.aliasName = 'type'
fm_type.outputField = type_name
fms.addFieldMap(fm_type)

fm_iucn = arcpy.FieldMap()
fm_iucn.addInputField(can, 'IUCN_CAT')
fm_iucn.addInputField(usa, 'IUCNcat')
iucn_name = fm_iucn.outputField
iucn_name.name = 'IUCN_cat'
iucn_name.aliasName = 'IUCN_cat'
fm_iucn.outputField = iucn_name
fms.addFieldMap(fm_iucn)

arcpy.Merge_management([can, usa], out_file, fms)


# join the original table to my working version using the uID, featureclasstofeatureclass
arcpy.env.qualifiedFieldNames = False
joined_table = arcpy.AddJoin_management('M08_mpa_20201124_FINALCOPY', 'uID_original', 'M05_mpa_MASTER_fieldsmapped', 'uID', 'KEEP_COMMON')
arcpy.CopyFeatures_management(joined_table, 'M09_mpa_joined')


## identify features that overlap
arcpy.CountOverlappingFeatures_analysis('M09_mpa_joined', 'OV01_mpa_overlap', 2, 'OV01_mpa_overlaptbl')
# (at this point, I noticed some slivers of things that weren't supposed to overlap. I fixed those manually, but they obvioulsy won't be reflected in previous versions)
arcpy.AddField_management('M09_mpa_joined', 'overlap', 'TEXT')
# go through overlap table, for each overlap, go through mpa fc and add a 'y' to the overlap field
with arcpy.da.SearchCursor('OV01_mpa_overlaptbl', ['ORIG_OID']) as cursor:
    for row in cursor:
        feature = row[0]
        with arcpy.da.UpdateCursor('M09_mpa_joined', ['OBJECTID', 'overlap']) as cursor:
            for row in cursor:
                if feature == row[0]:
                    row[1] = 'y'
                    cursor.updateRow(row)




# place a copy in the archive gdb


# for when its time to run simulations:
# identify ones that I most likely will not include in the final run (e.g. Olympic Coast and San Juans, since smaller more protected MPAs are embedded in them).
# For Canada, I could still run with overlapping polys, since I may want to know all of this connectivity. I think it is fine if the particle counts get doubled in that area. However, this would create some issues with settlement and assigning it to a certain polygon. Perhaps in the settlement version I will need one where the overlap areas are erased? I could do a UNION and then erase one of the overlap features.
