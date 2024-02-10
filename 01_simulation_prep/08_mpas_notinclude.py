# THERE ARE MANUAL STEPS in this script.
# BE SURE that datasets are backed up in mpas_archive.gdb before overwriting.

# Identify MPAs that should not be included in the analysis.
# 2021-06-08
# This is being done after the simulations. I somewhat recorded these decisions
# earlier, but I need this documented better. It will also be used to remove
# some MPAs that were included in the simulation but whose results should now
# be removed for any future analysis.

# This involves MANUALLY going through the original individual Canada and US
# protected areas and making notes on the ones that were not included.
# I want to end up with a list of protected areas that were:
#   (1) never included from the start because they are technically a terrestrial
#       protected area, or because they are in an area that is not resolved by
#       the hydrodynamic models,
#   (2) or, were included in the simulations, but should now be excluded from
#       analysis because they are in areas that are in areas that are not
#       resolved enough.

import arcpy
import os

base = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA'
gdb = os.path.join(base, 'mpas_notinclude.gdb')
arcpy.env.workspace = gdb

# NEP model extent
nep = os.path.join(base, 'mpas.gdb/A00_nep36_mb')
# High res coastline
coast_can = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_boundaries.gdb\Provinces'

# prep CARTS data
# I am using the original CARTS database because I need to include more than
# just the protected areas that are classified as Marine biome. Some terrestrial
# areas overlap the ocean. I did not inlcude these in my simulations, but I want
# to note where those areas are.
# I did not inlcude them because they likely aren't managed as MPAs and were not
# designed to protect marine resources.
# Clip to model extent
# Erase with coastline
# Dissolve by protected area to create multipart features
carts = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\1Data_BASE\BASE\BASE_protected_areas.gdb\CARTS_20171231'
arcpy.Clip_analysis(carts, nep, 'M000_mpa_can_1clip')
arcpy.Erase_analysis('M000_mpa_can_1clip', coast_can, 'M000_mpa_can_2erase')
arcpy.Dissolve_management(
    'M000_mpa_can_2erase', 
    'M001_mpa_can', 
    'PARENT_ID',
    [['NAME_E', 'FIRST'], ['BIOME', 'FIRST']])


# copy other datasets to gdb
# USA: M02_mpa_usa_erase is from my original MPA processing. It includes the
#      original protected areas clipped to the model extent and erased by the
#      high resolution coastline.
# m09 is the final version of my MPAs that I used in the Opendrift simulations.
mpa_usa = os.path.join(base, 'mpas.gdb/M02_mpa_usa_erase')
m09 = os.path.join(base, 'mpas.gdb/M09_mpa_joined')
arcpy.CopyFeatures_management(mpa_usa, 'M002_mpa_usa')
arcpy.CopyFeatures_management(m09, 'M003_M09')

# add attribute to each for noting exclusion
arcpy.AddField_management('M001_mpa_can', 'exclude', 'SHORT')
arcpy.AddField_management('M002_mpa_usa', 'exclude', 'SHORT')
arcpy.AddField_management('M003_M09', 'exclude', 'SHORT')

# MANUALLY:
# Go through each MPA in each dataset.
# For the original datasets, compare to the M09 dataset to see if it was
# included.
# For the M09 dataset, compare to the oceanographic model pts to see if it 
# should be excluded.
# Criteria:
#   Code as "1" for any that have no exposure to open ocean, are in narrow areas
#   and should definitely be excluded (ones that I've already drawn in M09 and
#   should be excluded).
#   Code as "2" if they are in very wide inlets and could maybe included.
#   Code as "3" if they are terrestrial and were not included.
#   Code as "4" if they are completely pelagic and were not included.


# Create a table of M09 with just the attribute and the uID
arcpy.TableToTable_conversion('M003_M09', gdb, 'M004_M09_table')
arcpy.DeleteField_management(
    'M004_M09_table',
    ['uID_original', 'partID', 'comments', 'OBJECTID_1', 'uID', 'name', 'provstate', 'gov_lev', 'mgmt', 'type', 'IUCN_cat', 'overlap', 'Shape_Length', 'Shape_Area']
    )

# Send a copy to the main mpa gdb and give it a name: M10_excludefromanalysis
main_mpa_gdb = os.path.join(base, 'mpas.gdb')
arcpy.Copy_management('M004_M09_table', os.path.join(main_mpa_gdb, 'M10_toexcludefromanalysis'))

# Send copies of the three datasets to mpa_archive.gdb
archive_gdb = os.path.join(base, 'mpas_archive.gdb')
arcpy.CopyFeatures_management('M001_mpa_can', os.path.join(archive_gdb, 'M001_mpa_can_TOEXCLUDE'))
arcpy.CopyFeatures_management('M002_mpa_usa', os.path.join(archive_gdb, 'M002_mpa_usa_TOEXCLUDE'))
arcpy.CopyFeatures_management('M003_M09', os.path.join(archive_gdb, 'M003_M09_TOEXCLUDE'))

# Then, create quick and dirty maps in a new arpx of the ones excluded.