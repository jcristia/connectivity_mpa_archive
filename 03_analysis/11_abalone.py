# Northern abalone case study
# Regenerate connectivity lines and destination rasters using just particles 
# that were released from high suitability habitat in MPAs and that settle on 
# moderate/high suitable habitat.

import arcpy
from arcpy.sa import *
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import pandas as pd
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import seaborn as sns

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity_chap2'
out_gdb = os.path.join(root, r'scripts\analysis_results\abalone.gdb')
arcpy.env.workspace = out_gdb

coastline = os.path.join(root, r'spatial\Coastline\coastline.gdb\landmask_FINAL')
mpas = 'M09_mpa_joined_M10exclude'
habitat = 'habitat_suitability_jamieson2004'


ncf = r'D:\MPA_connectivity_BACKUP\scripts\sim{}{}\outputs\nc\output_{}.nc'
origin_lat = os.path.join(root, r'cluster_results/scripts/sim{}{}/outputs/npy/lat_{}.npy')
origin_lon = os.path.join(root, r'cluster_results/scripts/sim{}{}/outputs/npy/lon_{}.npy')
dest_pts = os.path.join(root, r'cluster_results/scripts/sim{}{}/outputs/shp/dest_pts/dest_biology_pts_{}_pld{}.shp')

season = '05'
years = ['11', '14', '17']
plds = ['3', '7', '10']
mpa_groups = ['1', '2', '3']


#=================================
# Habitat suitability data

# Clip high habitat suitability polygons to MPAs
sel = arcpy.SelectLayerByAttribute_management(habitat, "NEW_SELECTION", "habitat_suitability = 'high'")
arcpy.Clip_analysis(sel, mpas, 'habitat_suitability_high_cliptoMPAs')

# Buffer all habitat polys similar to how I did it for the MPAs
# add unique ID
arcpy.CalculateField_management(habitat, 'uID', "!OBJECTID!", "PYTHON3", field_type="SHORT")
buff_fc = arcpy.Buffer_analysis(habitat, 'hs_TEMP1', '200 meters', dissolve_option='NONE')
# clip to coastline
clipped_fc = arcpy.Clip_analysis(buff_fc, coastline, 'hs_TEMP2')
# merge this with original
merged_fc = arcpy.Merge_management([habitat, clipped_fc], 'hs_TEMP3')
# dissolve by uID
diss_fc = arcpy.Dissolve_management(merged_fc, 'habitat_suitability_all_buff', 'uID')
arcpy.Delete_management(['hs_TEMP1', 'hs_TEMP2', 'hs_TEMP3'])

#=================================
# Get origin locations and ids for just the high suitability abalone polys 
# within MPAs.

for year in years:
    for group in mpa_groups:

        print(f"Getting origin locations for year {year} group {group}")

        dataset = nc.Dataset(ncf.format(year, season, group), "r")
        traj = dataset.variables["trajectory"]

        lats = np.load(origin_lat.format(year, season, group))
        lons = np.load(origin_lon.format(year, season, group))

        df = pd.DataFrame()
        df['lat'] = list(lats)
        df['lon'] = list(lons)
        df['traj_id'] = list(traj)
        dataset.close()
        del dataset

        x = np.array(np.rec.fromrecords(df.values))
        names = df.dtypes.index.tolist()
        x.dtype.names = tuple(names)
        sr = arcpy.SpatialReference(4326)
        arcpy.da.NumPyArrayToFeatureClass(x, os.path.join(arcpy.env.workspace, 'origin_pts_TEMP'), ('lon', 'lat'), sr)

        arcpy.Project_management('origin_pts_TEMP', 'origin_pts_PROJECT', arcpy.SpatialReference(3005))
        arcpy.Clip_analysis('origin_pts_PROJECT', 'habitat_suitability_high_cliptoMPAs', 'origin_pts_{}{}_{}'.format(year, season, group))
        arcpy.Delete_management('origin_pts_TEMP')
        arcpy.Delete_management('origin_pts_PROJECT')

#=================================
# Get destination points that match the IDs of the origin points
# Get just the ones that settle on suitable habitat (moderate and highly 
# suitable habitat). Use a buffered version of the habitat polygons.

arcpy.env.qualifiedFieldNames = False

for year in years:
    for group in mpa_groups:
        for pld in plds:

            print(f"Destpts for year {year} group {group} pld {pld}")

            dpts = dest_pts.format(year, season, group, pld)
            opts = f'origin_pts_{year}{season}_{group}'

            arcpy.MakeFeatureLayer_management(dpts, 'temp_lyr')
            joined_fc = arcpy.management.AddJoin('temp_lyr', 'traj_id', opts, 'traj_id', 'KEEP_COMMON')
            arcpy.CopyFeatures_management('temp_lyr', 'temp_fc')

            outfc = f'destpts_{year}{season}_{group}_pld{pld}'
            arcpy.Clip_analysis('temp_fc', 'habitat_suitability_all_buff', outfc)

            arcpy.Delete_management(['temp_lyr', 'temp_fc'])

#=================================
# Create connectivity lines between MPAs based on habitat connectivity

sg = pd.DataFrame.spatial.from_featureclass(mpas)

for year in years:
    for group in mpa_groups:
        for pld in plds:

            print(f"Conn lines for year {year} group {group} pld {pld}")

            time_step_output = 0.5
            pld = int(pld)
            pld_int = int((pld * 24) / time_step_output)

            # read in destination points as df
            dpts = f'destpts_{year}{season}_{group}_pld{pld}'
            field_names = [i.name for i in arcpy.ListFields(dpts) if i.type != 'OID']
            cursor = arcpy.da.SearchCursor(dpts, field_names)
            od = pd.DataFrame(data=[row for row in cursor], columns=field_names)

            # on od, select particles where time_int_s minus time_int is less than or equal to PLD
            od_pld = od[(od.time_int - od.time_int_s <= pld_int)]

            # get each unique combination of originID and destID and get count of particles that survived
            od_unique = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID_part','dest_id']).size().reset_index(name='Freq')
            # how to read this:
            # first we select from od the ones that settled and survived
            # then we groupby unique combinations of uID and dest_id
            # then we get the count of those unique combinations
            # this normally makes uID the index and doesn't have a column name for count 
            # (the series we created), so we reset index and give the count a column name

            # df of time interval where first settlement occurred
            df_time_int = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID_part','dest_id'])['time_int'].min().reset_index(name='time_int')

            # set up for creating self connection lines. Size of circle lines based on amount settled and average area of all patches.
            def CircleCoords(xLeft, yCenter, r, n): # credit: MGET. Also, see circle_coords.xlsx for explanation of equation.
                return [(xLeft + r - math.cos(2*math.pi/n*x)*r, math.sin(2*math.pi/n*x)*r + yCenter) for x in range(n+1)]
            # min and max quantities used for normalization
            quantity_min = od_unique[od_unique.dest_id == od_unique.uID_part].Freq.min()
            quantity_max = od_unique[od_unique.dest_id == od_unique.uID_part].Freq.max()
            # get average area
            area_mean = sg.spatial.area / len(sg)
            # get radius of a circle with this area
            radius = math.sqrt(area_mean/math.pi)

            # for each unique combinaton create line from centroids
            connection_lines = pd.DataFrame(columns=['from_id','to_id','quantity','totalori','prob','time_int', 'line', 'pld'])
            conn_i = 0
            for row in od_unique.itertuples(index=False):
                # get total amount of particles released from patch
                total = od.uID_part[od.uID_part ==  row[0]].value_counts().values[0]

                # time interval where first settlement occurred
                time_int = df_time_int[(df_time_int.uID_part == row[0]) & (df_time_int.dest_id == row[1])]['time_int'].values[0]

                # get centroid of from and to patches
                centroid_origin = sg[sg.uID_20201124 == row[0]].spatial.centroid
                centroid_dest = sg[sg.uID_20201124 == row[1]].spatial.centroid

                if row[0] != row[1]:
                    geom_line = [centroid_origin, centroid_dest]
                else:
                    # normalize the quantites to 0.5 - 1 range (or I can do 0-1 but then the smallest one won't show up)
                    #quantity_norm = 0.5 * (row[2] - quantity_min) / float(quantity_max - quantity_min) + 0.5
                    quantity_norm = (row[2] - quantity_min) / float(quantity_max - quantity_min)
                    if quantity_norm == 0 or np.isnan(quantity_norm):
                        quantity_norm = 0.05
                    radius_adj = radius * quantity_norm
                    geom_line = CircleCoords(centroid_origin[0], centroid_origin[1], radius_adj, 90)

                connection_lines.loc[conn_i] = [row[0],row[1],float(row[2]),float(total),row[2]/float(total), time_int,geom_line, pld]
                conn_i += 1

            connection_lines['date_start'] = f'{year}{season}'

            #create blank fc, use existing lines fc as template to get fields
            template = os.path.join(root, 'cluster_results/scripts/COMBINED.gdb/connectivity_1105_pld3')
            out_fc = f'connectivity_{year}{season}_{group}_pld{pld}'
            arcpy.CreateFeatureclass_management(
                arcpy.env.workspace,
                out_fc,
                'POLYLINE',
                template=template,
                spatial_reference=template
                )

            field_names = ['from_id','to_id','quantity','totalori','prob','time_int', 'pld', 'date_start', 'SHAPE@']
            with arcpy.da.InsertCursor(out_fc, field_names) as cursor:

                for row in connection_lines.itertuples(index=False):
                    polyline = arcpy.Polyline(arcpy.Array([arcpy.Point(*coords) for coords in row[6]]))
                    cursor.insertRow([row[0], row[1], row[2], row[3], row[4], row[5], row[7], row[8], polyline])



#=================================
# Combine connectivity lines by PLD
# Average over time periods

for pld in plds:
    fcs = arcpy.ListFeatureClasses("connectivity*pld{}".format(pld))
    arcpy.Merge_management(fcs, f'connectivity_ALL_pld{pld}')

def mean_cust_denom(x):
    s = x.sum()
    m = s/float(len(plds))
    return m

for pld in plds:

    fc = f'connectivity_ALL_pld{pld}'
    sedf = pd.DataFrame.spatial.from_featureclass(fc)

    # average
    df_group = sedf.groupby(['from_id', 'to_id']).agg(
        freq = ('from_id', 'count'),
        prob_avg = ('prob', mean_cust_denom),
        totalori = ('totalori', 'sum'),
        totquant = ('quantity', 'sum'),
        date_start = ('date_start', 'first'),
        # SHAPE = ('SHAPE', 'first'), # doesn't work if you include shape field
        pld = ('pld', 'first')
        )
    df_group = df_group.reset_index()

    # join back to original to get geometry
    sedf_l = sedf[['from_id', 'to_id', 'SHAPE']]
    df_merge = df_group.merge(sedf_l, on=['from_id', 'to_id'])
    df_nodupes = df_merge.drop_duplicates(subset=['from_id', 'to_id'], keep='first')

    #output
    df_nodupes.spatial.to_featureclass(location=os.path.join(arcpy.env.workspace, f'connectivity_ALLavg_pld{pld}'))
    

#=================================
# Combine destination particles by PLD

for pld in plds:
    fcs = arcpy.ListFeatureClasses("destpts*pld{}".format(pld))
    arcpy.Merge_management(fcs, f'destpts_ALL_pld{pld}')


#=================================
# Raster analysis to get particles that settle on unprotected coastline
# Do just for moderate/high suitable habitat destination locations.

# Delete Identical points
# Where MPAs overlap, points get duplicated, which is what I needed to create the
# correct connection lines, but I don't want to double count them here.
for pld in plds:
    print(f'deleting identical points for pld{pld}')
    arcpy.CopyFeatures_management(f'destpts_ALL_pld{pld}', f'destpts_delidentical_pld{pld}')
    arcpy.DeleteIdentical_management(f'destpts_delidentical_pld{pld}', ['traj_id', 'uID_part', 'date_start'])

# Filter for particles that are not self connections and that survive until settlement.
for pld in plds:
    print(f'selecting by attribute for pld {pld}')
    arcpy.MakeFeatureLayer_management(f'destpts_delidentical_pld{pld}', 'temp_lyr')
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'NEW_SELECTION', '"uID_part" <> "dest_id"')
    # so the way I have dest_pts now, the max that time_int can be is the 
    # time_int_s + pld_int. This would mean that the particle was still drifting
    # at the end of that PLD.
    # So now, remove where mortality happened before stranding.
    # NOTE ON PRECOMPETENCY PERIOD and ABALONE
    # I'm allowing ones to live if they strand BEFORE mortality. This is how I 
    # did it with the generic MPA runs, and it is also how I calculate the 
    # connectivity lines above.
    # However, I recognize that for abalone settelment may not occur until day 3
    # and therefore, larval mortality could still be occurring.
    # To not complicate things, I'm not going to keep applying mortality if they
    # strand early though. From what I read in the abalone literature I cite, 
    # if they encounter suitable habitat then they may move benthically and 
    # somewhat maintain their position until they are developed and able to
    # securely settle. Therefore, it's reasonable to believe that mortality 
    # would be much lower if they strand on suitable habitat.
    arcpy.SelectLayerByAttribute_management('temp_lyr', 'REMOVE_FROM_SELECTION', '"mortstep" < "time_int" and "mortstep" <> -1')
    arcpy.CopyFeatures_management('temp_lyr', f'destpts_selattr_pld{pld}')
    arcpy.Delete_management('temp_lyr')


# UGH, so because I drew my habitat polygons out to 30m, and some MPAs only 
# extend to 15m depth, I'm left with a lot of particles on the seaward side of 
# MPAs. In reality, these would effectively be "protected" particles and not
# something I should count as recruiting to unprotected areas.
# I don't want to edit the habitat polys, since that would be inconsistent, so
# manually create polys to mask these particles out.
# Create empty fc
# Manually trace polys that are in front of MPAs.
# Select points by location, inverse selection and those will be my points to
# rasterize.
arcpy.CreateFeatureclass_management(arcpy.env.workspace, 'destpts_mask', 'POLYGON')
# Unfortunately, if I am considering these protected, then they should be
# considered part of the connection strengths, but its a very minimal amount.
for pld in plds:
    arcpy.MakeFeatureLayer_management(f'destpts_selattr_pld{pld}', 'templyr')
    arcpy.SelectLayerByLocation_management('templyr', 'INTERSECT', 'destpts_mask', invert_spatial_relationship='INVERT')
    arcpy.CopyFeatures_management('templyr', f'destpts_selloc_pld{pld}')
    arcpy.Delete_management('templyr')

# convert destination particles to raster
snapras = os.path.join(root, r'spatial\Coastline\coastline.gdb\landmask_NN')
arcpy.env.snapRaster = snapras
# snap to the landmask raster which is 500m cell size. However, only snap to it the first time
# and then change the snap raster to the first one created. I am creating rasters with 1000m
# cell size, so if I keep snapping to the landmask then they don't always get created on
# the same grid
newSnapRas = True
for pld in plds:
    fc = f'destpts_selloc_pld{pld}'
    outname = f'destpts_recruitCount_pld{pld}'
    arcpy.PointToRaster_conversion(fc, '', outname, 'COUNT', cellsize=1000)
    if newSnapRas:
        arcpy.env.snapRaster = outname
        newSnapRas = False


# I want to know recruitment to only unprotected areas, so...
# remove cells that overlap the BUFFERED MPA polygons (turn to null values)

# rasterize mpa buff
mpa_buff = os.path.join(root, r'spatial\MPA\mpas_shp_release\mpas.gdb\MB06_mpa_buff_FINAL')
mpa_buff_out = os.path.basename(mpa_buff)
arcpy.CopyFeatures_management(mpa_buff, mpa_buff_out)
arcpy.AddField_management(mpa_buff_out, 'priority', 'SHORT')
with arcpy.da.UpdateCursor(mpa_buff_out, ['priority']) as cursor:
    for row in cursor:
        row[0]=1
        cursor.updateRow(row)
snapras = 'destpts_recruitCount_pld{}'.format(plds[0])
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
mpas_buff_isnull = mpa_buff_out + '_ras_isnull'
for pld in plds:
    ras = f'destpts_recruitCount_pld{pld}'
    outras = SetNull(mpas_buff_isnull, ras, 'Value = 0')
    outras.save(f'destpts_recruitCount_rmMPAs_pld{pld}')


# Quantifying unprotected habitat

# Unlike the coastwide analysis that was not habitat based and where I only
# quantified settlement right on the coastline (rasterized just the line
# feature), here, I want to inlcude all unprotected cells down to 30m (except 
# those noted above).

# Rasterize all unprotected habitat (except for features noted above to exclude)

# erase habitat by MPAs
arcpy.Erase_analysis(habitat, mpas, 'temp_haberasempas')
# erase habitat by destpt mask polygon
arcpy.Erase_analysis('temp_haberasempas', 'destpts_mask', 'temp_haberasemask')
# Annoyingly there a few slivers, which results in creating raster cells over
# protected areas. Delete these by uID.
uIDs_delete = [4,13,14,24,30,31,58]
with arcpy.da.UpdateCursor('temp_haberasemask', ['uID']) as cursor:
    for row in cursor:
        if row[0] in uIDs_delete:
            cursor.deleteRow()
# rasterize like how I did for MPA buff, but I need to exclude any overlapping
# MPA buff.
arcpy.AddField_management('temp_haberasemask', 'priority', 'SHORT')
with arcpy.da.UpdateCursor('temp_haberasemask', ['priority']) as cursor:
    for row in cursor:
        row[0]=1
        cursor.updateRow(row)
snapras = 'destpts_recruitCount_pld{}'.format(plds[0])
arcpy.env.snapRaster = snapras
arcpy.env.cellSize = snapras
arcpy.PolygonToRaster_conversion('temp_haberasemask', 'OBJECTID', 'temp_habrast', 'MAXIMUM_AREA', priority_field='priority')

# This overlaps with the MPA raster. Therefore, set the cells in the habitat 
# raster to null where it overlaps.
outras = SetNull(mpas_buff_isnull, 'temp_habrast', 'Value = 0')
outras.save(f'habitat_suitability_rast_nompasmask')

arcpy.Delete_management(['temp_habrast', 'temp_haberasemask', 'temp_haberasempas'])


#=================================
# Quantify recruitment

# total cell count
total_cells = 0
with arcpy.da.SearchCursor('habitat_suitability_rast_nompasmask', ['Value', 'Count']) as cursor:
    for row in cursor:
        total_cells += row[1]

# get coastline cells covered for each pld and threshold level
thresholds = list(range(0, 1010, 10))
df = pd.DataFrame()
for pld in plds:
    print(f'processing pld {pld}')
    fc = f'destpts_recruitCount_rmMPAs_pld{pld}'
    for thresh in thresholds:
        # set the recruit raster to 0 where less than threshold and 1 where above
        outCon = Con(fc, 0, 1, f"Value < {thresh}")
        # save the first one for reference:
        if thresh == 0:
            outCon.save(f'destpts_pld{pld}_thresh{thresh}')
        outCon.save('temp_ras')
        # count cells and divide by total
        with arcpy.da.SearchCursor('temp_ras', ['Value', 'Count']) as cursor:
            for row in cursor:
                if row[0]==1:
                    cell_count = row[1]
        arcpy.Delete_management('temp_ras')

        ### Total quantity settled by threshold ###
        # set the recruit raster to 0 where less than the threshold
        outCon2 = Con(fc, 0, fc, f'Value < {thresh}')
        # using the attribute table, multiply the count by the value and add up
        outCon2.save(f'temp_ras2')
        total_particles = 0
        with arcpy.da.SearchCursor('temp_ras2', ['Value', 'Count']) as cursor:
            for row in cursor:
                total_particles += row[0] * row[1]
        arcpy.Delete_management('temp_ras2')

        # Get total quantity released
        fc_2 = f'destpts_delidentical_pld{pld}'
        total_released = arcpy.GetCount_management(fc_2)[0]

        # append to dataframe
        df = df.append(
            {'pld':pld, 
            'threshold':str(thresh), 
            'percent_cover':(cell_count/total_cells*100),
            'total_recruits': total_particles},
            ignore_index=True
            )

# Get total amount released
total_rel = 0
for year in years:
    for group in mpa_groups:
        fc = f'origin_pts_{year}{season}_{group}'
        total = arcpy.GetCount_management(fc)[0]
        total_rel += int(total)

df['total_released'] = total_rel
df['percent_settled'] = df.total_recruits/df.total_released*100

df.to_csv(os.path.join(root, 'scripts/analysis_results', 'recruit_totals_abalone.csv'))

#=================================
# Plots
df = pd.read_csv(os.path.join(root, 'scripts/analysis_results', 'recruit_totals_abalone.csv'))

df['PLD'] = df.pld.astype('int')
df['PLD'] = df.PLD.astype('str') # seaborn hue needs to be a string
df = df[df.threshold <= 200]

sns.set(rc = {'figure.figsize':(17,8)})
sns.set_style('white')
#sns.set_context('paper', font_scale=1.5, rc={"lines.linewidth": 2})
sns.set_context('paper', font_scale=2, rc={"lines.linewidth": 3.5})
fig,axs = plt.subplots(ncols=2)
f = sns.lineplot(
    data = df,
    x = 'threshold',
    y = 'percent_cover',
    hue = 'PLD',
    ax=axs[0]
)
axs[0].legend('', frameon=False)
f.set(xlabel='Threshold (particles settled per raster cell)', ylabel='% of unprotected suitable abalone habitat\nreceiving particles from protected habitat')
g = sns.lineplot(
    data = df,
    x = 'threshold',
    y = 'percent_settled',
    hue = 'PLD',
    ax=axs[1]
)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='PLD (days)')
g.set(xlabel='Threshold (particles settled per raster cell)', ylabel='% of total particles released that settled on\nunprotected suitable abalone habitat')
f.figure.savefig('fig_abalone_recruits_settled.svg')

# PLD of 3 is higher for % of unprotected habitat because a lot of habitat
# patches were directly adjacent to protected ones, and therefore longer PLDs
# would drift over these patches.
# For where these cross over with increasing threshold, my guess is that many
# for PLD 3 do drift out of any patches and don't have enough time to get back.


#=================================
# Additional metric to state in methods/results:
# % of high suitability habitat in region covered by MPAs

sedf_all = pd.DataFrame.spatial.from_featureclass(habitat)
all_high_area = sedf_all.SHAPE.geom.area[sedf_all.habitat_suitability=='high'].sum()
sedf_clip = pd.DataFrame.spatial.from_featureclass('habitat_suitability_high_cliptoMPAs')
mpa_high_area = sedf_clip.SHAPE.geom.area.sum()
percent = mpa_high_area / all_high_area * 100
# percent: 49.02%

