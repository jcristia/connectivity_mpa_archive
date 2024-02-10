# (1) Calculate % of coast that receives recruits
# (2) Calculate number of recruits that settle on coast
# Break this out by PLD and a Threshold level
# Only do for the BC coast and only for the "outer" coast.
# Between these to plots, I can get a measurement of relative coverage and 
# overall recruitment.
# For the % of coast, there is a small amount of error due to some particle
# ending positions being a pixel too far inland, or just off the coast. I did
# a visual survey of this and the effect should be small.


import arcpy
from arcpy.sa import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity_chap2'
out_gdb = os.path.join(root, r'scripts\analysis_results\raster_recruit_summary.gdb')
recruit_rasters_gdb = os.path.join(root, r'cluster_results\scripts\DEST_RAST.gdb')
coastline = os.path.join(root, r'spatial\Coastline\coastline.gdb\landmask_FINAL')
arcpy.env.workspace = out_gdb


### rasterize coastline ###
arcpy.PolygonToLine_management(coastline, 'coastline_line')
# clip to poly that excludes areas of the coast not resolved as well as the US 
# coast. This was drawn manually.
arcpy.Clip_analysis('coastline_line', 'coast_clip_poly', 'coastline_line_clip')
# Dummy value for raster
arcpy.AddField_management('coastline_line_clip', 'dummy', 'SHORT')
arcpy.CalculateField_management('coastline_line_clip', 'dummy', 1)
snapras = os.path.join(recruit_rasters_gdb, 'recruit_count_ALL_1')
arcpy.env.snapRaster = snapras
arcpy.env.cellSize = snapras # use the same cell size otherwise you might miss overlap and/or double count
arcpy.PolylineToRaster_conversion('coastline_line_clip', 'dummy', 'coastline_rast')


### % of coastline receiving recruits ###

# total cell count. They are all 1s.
with arcpy.da.SearchCursor('coastline_rast', ['Value', 'Count']) as cursor:
    for row in cursor:
        if row[0]==1:
            total_cells = row[1]

# get coastline cells covered for each pld and threshold level
plds = ['1', '3', '7', '10', '21', '30', '40', '60']
thresholds = list(range(0, 1010, 10))
df = pd.DataFrame(columns={'pld', 'threshold', 'perc_cover'})
df = pd.DataFrame()
for pld in plds:
    print(f'processing pld {pld}')
    fc = os.path.join(recruit_rasters_gdb, f'recruit_count_ALL_{pld}')
    for thresh in thresholds:
        # set the recruit raster to 0 where less than threshold and 1 where above
        outCon = Con(fc, 0, 1, f"Value < {thresh}")
        # set the recruit raster to null where the coastline raster is null
        outras = SetNull(outCon, 'coastline_rast', 'Value = 0')
        # save the first one for reference:
        if thresh == 0:
            outras.save(f'recruit_setnull_{pld}_{thresh}')
        outras.save('temp_ras')
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

        # append to dataframe
        df = df.append(
            {'pld':pld, 
            'threshold':str(thresh), 
            'perc_cover':(cell_count/total_cells),
            'total_recruits': total_particles}, 
            ignore_index=True
            )

df.to_csv(os.path.join(root, 'scripts/analysis_results', 'recruit_totals.csv'))
df = pd.read_csv(os.path.join(root, 'scripts/analysis_results', 'recruit_totals.csv'))


#### plots ####

df['PLD'] = df.pld.astype('int')
df['PLD'] = df.PLD.astype('str') # seaborn hue needs to be a string

df = df[df.threshold <= 200]
df['percent_cover'] = df.perc_cover *100

sns.set()
sns.set_style('white')
sns.set_context('paper', font_scale=1.5, rc={"lines.linewidth": 2})
f = sns.lineplot(
    data = df,
    x = 'threshold',
    y = 'percent_cover',
    hue = 'PLD',
    #palette = 'husl'
)
plt.legend([],[], frameon=False)
f.set(xlabel='Threshold (particles settled per raster cell)', ylabel='% of coast receiving particles')
f.figure.savefig('recruits_perc_cover.jpg')


# get total particle count to compare next plot to
mpas = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity_chap2\spatial\MPA\mpas_shp_release\mpa_.shp'
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity_chap2\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'
field_names = [i.name for i in arcpy.ListFields(mpas) if i.type != 'OID']
cursor = arcpy.da.SearchCursor(mpas, field_names)
df_mpas = pd.DataFrame(data=[row for row in cursor], columns=field_names)
field_names = [i.name for i in arcpy.ListFields(mpas_exclude) if i.type != 'OID']
cursor = arcpy.da.SearchCursor(mpas_exclude, field_names)
df_mpas_ex = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_mpas = df_mpas.merge(df_mpas_ex, left_on='uID_202011', right_on='uID_20201124')
df_mpas = df_mpas[df_mpas.exclude != 1]
df_mpas = df_mpas[(df_mpas.govID.str.startswith('5')) | (df_mpas.govID.str.startswith('7'))] # exclude US mpas
total_part = df_mpas.part_num.sum()
total_part = total_part * 84 # releases
total_part = total_part * 9 # time periods

# interpretation:
# there were ~17 million particles released
# compared the total amount settled in figure 2 to this number
# It is a rough total estimate though. Particles from excluded MPAs and MPAs in
# the USA could have drifted into ok territory.

df['percent_settled'] = df.total_recruits / total_part * 100

f = sns.lineplot(
    data = df,
    x = 'threshold',
    y = 'percent_settled',
    hue = 'PLD',
    #palette = 'husl'
)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='PD')
f.set(xlabel='Threshold (particles settled per raster cell)', ylabel='% of total particles release that settled on coast')
f.figure.savefig('recruits_total_settled.jpg')
# if legend gets cut off, just save it manually from visual studio



# PLOT THEM ALL TOGETHER!!!
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
f.set(xlabel='Threshold (particles settled per raster cell)', ylabel='% of coast receiving particles')
g = sns.lineplot(
    data = df,
    x = 'threshold',
    y = 'percent_settled',
    hue = 'PLD',
    ax=axs[1]
)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='PLD (days)')
g.set(xlabel='Threshold (particles settled per raster cell)', ylabel='% of total particles released that settled on coast')
f.figure.savefig('fig_06ab_recruits_settled.svg')



# Additional analysis following CJFAS review
# Calculate the percent of coastline that is not protected to add context to the
# results.

# total cell count. They are all 1s.
with arcpy.da.SearchCursor('coastline_rast', ['Value', 'Count']) as cursor:
    for row in cursor:
        if row[0]==1:
            total_cells = row[1]
# total cells: 10025

mpa_rast = os.path.join(recruit_rasters_gdb, 'Con_MB06_mpa1')

# MPA rast is already 0
# Coast rast is 1
outras = Raster('coastline_rast') * Raster(mpa_rast)
outras.save('coastline_mpas')

# count cells
with arcpy.da.SearchCursor('coastline_mpas', ['Value', 'Count']) as cursor:
    for row in cursor:
        if row[0]==0:
            cell_count = row[1]
# cell count: 3384

# Proportion of coast that is not protected:
coast_not_protected = (total_cells - cell_count) / total_cells
# coast_not_protected: 66%
