# script to format the mpa shapefiles that are used for particle release in 
# Opendrift

# env: arcgispro-py3-clone_MPACONNECTIVITY

import arcpy

mpas_fc = (
    r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA'
    '\mpas.gdb\M08_mpa_20201124_FINALCOPY')
mpas_shp_folder = (
    r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA'
    '\mpas_shp_release')
gdb = 'mpas.gdb'
fc = os.path.basename(os.path.normpath(mpas_fc))
arcpy.env.workspace = os.path.join(mpas_shp_folder, gdb)



if not arcpy.Exists(arcpy.env.workspace):
    arcpy.CreateFileGDB_management(mpas_shp_folder, gdb)
if not arcpy.Exists(fc):
    arcpy.CopyFeatures_management(mpas_fc, fc)
    arcpy.AddField_management(fc, 'part_num', 'LONG')

############################################
# equation for calculating particles to release

# get area min and max
with arcpy.da.SearchCursor(fc, ['Shape_Area']) as cursor:
    area_max = max(cursor)
    cursor.reset()
    area_min = min(cursor)

# min: 9
# max: 1,224,332,320 (Olympic Coast)
# 2nd largest: 355,986,778
# mean: 15,814,497
# median: 722,623
# the largest Canadian MPA is Gwaii Haanas (well, one part of it): 310,867,762

# similar to my seagrass chapter, I will have 84 releases (release every 4 hours
# for two weeks)
# I want my base amount of particles per release to be 10
# If I do things propotionally, that equals 1.3 billion per release for the
# largest MPA. Obviously way too much.

####### Notes for how I calculated particle amount to release:

"""
I can't do a linear relationship between particle count and area. There is a 
massive range in order of magnitude of patch sizes. I wanted to make sure I 
released enough for the smaller patches without having to release a massive 
amount for larger patches.

I explored some simple clean equations, but none of them worked that well 
(log, ln, square root, etc). To get the relationship that I want, I need to 
define a few points. I the end this will be the easiest to justify.

I set the amounts I want to release for the min, max, and median patch size. 
For min, I wanted to release at least 1 particle per release. This would be 84
total, which means my min connection strength would be 0.012. For max, I 
examined the few largest patches. They are massive. I set this at 1,000 (I can 
then detect probabilities as low as 0.001). This would be for Olympic 
Coast, which is an order of magnitude larger than the others. I'm simply not 
as concerned with this one anyways. For the median, which is many order of 
magnitude less than the max, I set this at 50 particles per release. Based 
on a visual assessment of these patches, I think this is adequate.

I put these 3 data points into mycurvefit.com (quick and dirty), and fit a power
 curve. This allows for quick growth at the beginning and then growth slows.

When compared to Masson 2015, I am releasing than them per 
unit area (except for the largest patch). They also didn't offer any 
justification for their particle count. I can say something like "it exceeds the
particles released from past studies on the BC coast (e.g. ...; ...)"
Masson: 0.0000022 particles per meter squared
Mine:
smallest   0.11/1 m2
median 0.000069/1 m2 (if I release 50)
largest 0.00000082/1 m2

However, I will eventually do sensitivity tests on the smallest, median, 
and largest patches (or patches near these sizes that are spaced evenly 
throughout). See my notes in Evernote about particle count sensitivity tests.

equation from mycurvefit.com (yes, quick and dirty)
y = 0.2180685x^0.4028882
"""

############################################

# calculate particles based on area
with arcpy.da.UpdateCursor(fc, ['Shape_Area', 'part_num']) as cursor:
    for row in cursor:
        particles = 0.2180685 * (row[0] ** 0.4028882)
        row[1] = particles
        cursor.updateRow(row)

############################################

# split out mpas into individual shapefiles
# (NO LONGER DOING THIS, but I will keep this code for now)

# # export to shapefile
# shp = os.path.join(mpas_shp_folder, 'mpa_.shp')
# arcpy.CopyFeatures_management(fc, shp)

# # create folder if it doesn't exists
# mpas_shp = os.path.join(mpas_shp_folder, 'mpas_shp')
# if not os.path.exists(mpas_shp):
#     os.makedirs(mpas_shp)

# # split
# arcpy.SplitByAttributes_analysis(shp, mpas_shp, 'uID_202011')


#############################################

# grouping
# for Opendrift, I tired simulating 1 feature at a time, thinking that this
# would be faster, but it was actually quite slow. It is better to simulate many
# patches at once. Therefore, I will manually separate them into nearby groups
# (assuming they will run faster if not as much of the reader is called),
# but keep it as one shapefile.

# In my seagrass chapter I did 4,000 total particles per release. Stick with
# that. Add up polygons and get about 4,000 particles.

# add attribute
arcpy.AddField_management(fc, 'grouping', 'SHORT')

# NOW:
# do groupings visually
# do it on a copied version:

#arcpy.CopyFeatures_management(fc, 'M08_mpa_20210120_groupings')

# then manually export to shape (name: mpa_.shp)